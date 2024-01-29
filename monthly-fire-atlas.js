/**
 * This is an example code of the following paper:
 * 
 * Alana K. Neves, Manuel L. Campagnolo, João M.N. Silva, José M.C. Pereira, A Landsat-based atlas of monthly burned area for Portugal, 1984–2021, International Journal of Applied Earth Observation and Geoinformation, Volume 119, 2023, 103321, ISSN 1569-8432, https://doi.org/10.1016/j.jag.2023.103321.
 * 
 * Since the processing is performed separately for each year, 
 * the code below is only for the year 2016.
 * 
*/

//Import mainland Portugal shapefile as feature collection
var PT_int = ee.FeatureCollection('users/alanakneves/Portugal-limit-continente');
//Collect the feature collection geometry
var geom = PT_int.geometry();

print(geom);

Map.centerObject(geom, 6);

//Start and end dates of analysis
var start = ee.Date('2016-01-01');
var end = ee.Date('2016-12-31');

// FILTER LANDSAT 7 AND 8 IMAGES BY DATE, AREA OF INTEREST AND CLOUD COVER

var filteredCollection_all_8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filter(ee.Filter.bounds(geom))
  .filterMetadata('CLOUD_COVER', 'less_than', 40)
  .filterDate(start, end).sort('system:time_start');
  
var filteredCollection_all_7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filter(ee.Filter.bounds(geom))
  .filterMetadata('CLOUD_COVER', 'less_than', 40)
  .filterDate(start, end).sort('system:time_start');
  
print('Filtered collection L8  - all year: ', filteredCollection_all_8);
print('Filtered collection L7  - all year: ', filteredCollection_all_7);

// APPLY FOCAL MEAN TO FILL L7 IMAGE GAPS

var filtered_focal_mean = function(img) {
  var img_fill = img.focalMean(1, 'square', 'pixels', 8).blend(img);
  return img_fill.int().copyProperties(img,['system:time_start']);
};

var l7_filtered = filteredCollection_all_7.map(filtered_focal_mean);

print('L7 after focal mean', l7_filtered);

// CHECK IF L7 GAP FILLING WORKED

var list_l7 = filteredCollection_all_7.sort('system:time_start').toList(filteredCollection_all_7.size());
print('List of all L7 images', list_l7);
var l7_single = ee.Image(list_l7.get(51));
Map.addLayer(l7_single, {bands: ['SR_B7', 'SR_B4', 'SR_B3'], min: 3000, max: 20000}, 'L7 before gap filling',false);

var list_l7_filt = l7_filtered.sort('system:time_start').toList(l7_filtered.size());
print('List of all L7 images after focal mean', list_l7_filt);
var l7_single_ = ee.Image(list_l7_filt.get(51));
Map.addLayer(l7_single_, {bands: ['SR_B7', 'SR_B4', 'SR_B3'], min: 3000, max: 20000}, 'L7 after gap filling',false);

//CREATE A DOY BAND FOR EACH IMAGE IN COLLECTION

var addDate = function(image){
  var doy = image.date().getRelative('day', 'year');
  var doyBand = ee.Image.constant(doy).uint16().add(1).rename('doy');
  doyBand = doyBand.updateMask(image.select('SR_B5').mask());
  
  return image.addBands(doyBand);
};

// Incorporate DOY bands into the original collections

var withDate_l7 = l7_filtered.map(addDate).sort('system:time_start');
print('L7 Collection with doy band: ', withDate_l7);

var withDate_l8 = filteredCollection_all_8.map(addDate).sort('system:time_start');
print('L8 Collection with doy band: ', withDate_l8);

//RENAME L7 AND L8 BANDS TO MERGE THEM IN AN UNIQUE COLLECTION

////L7
var rename_l7 = function(img) {
    var select = img.select('SR_B3','SR_B4','SR_B7','ST_B6','doy','QA_PIXEL','QA_RADSAT')
                    .rename(['Red','NIR','SWIR2','Thermal','doy','QA_PIXEL','QA_RADSAT'])
                  .set('system:time_start', img.get('system:time_start'));
    return select;
};

var l7_renamed = withDate_l7.map(rename_l7);
print('L7 with renamed bands:', l7_renamed);

var list_l7_ren = l7_renamed.sort('system:time_start').toList(l7_renamed.size());
var l7_ren_test = ee.Image(list_l7_ren.get(51));
Map.addLayer(l7_ren_test, {bands: ['SWIR2', 'NIR', 'Red'], min: 3000, max: 20000}, 'L7 with renamed bands',false);

////L8
var rename_l8 = function(img) {
    var select = img.select('SR_B4','SR_B5','SR_B7','ST_B10','doy','QA_PIXEL','QA_RADSAT')
                    .rename(['Red','NIR','SWIR2','Thermal','doy','QA_PIXEL','QA_RADSAT'])
                  .set('system:time_start', img.get('system:time_start'));
    return select;
};

var l8_renamed = withDate_l8.map(rename_l8);
print('L8 with renamed bands:', l8_renamed);

var list_l8_ren = l8_renamed.sort('system:time_start').toList(l8_renamed.size());
var l8_ren_test = ee.Image(list_l8_ren.get(10));
Map.addLayer(l8_ren_test, {bands: ['SWIR2', 'NIR', 'Red'], min: 3000, max: 20000}, 'L8 with renamed bands',false);

// MERGE L7 AND L8 COLLECTIONS

var filteredCollection_all = l8_renamed.merge(l7_renamed).sort('system:time_start');
print('Filtered collection - L7 and L8 merged: ', filteredCollection_all);

// CREATE MOSAICS WITH IMAGES FROM THE SAME DATE

function mosaicByDate(imcol){
  // imcol: An image collection
  // returns: An image collection
  var imlist = imcol.toList(imcol.size());

  var unique_dates = imlist.map(function(im){
    return ee.Image(im).date().format("YYYY-MM-dd");
  }).distinct();

  var mosaic_imlist = unique_dates.map(function(d){
    d = ee.Date(d);

    var im = imcol
      .filterDate(d, d.advance(1, "day"))
      .mosaic();

    return im.set(
        "system:time_start", d.millis(), 
        "system:id", d.format("YYYY-MM-dd"));
  });

  return ee.ImageCollection(mosaic_imlist);
}

var col_mosaic = mosaicByDate(filteredCollection_all).sort('system:time_start');

print('Mosaics of the same date', col_mosaic);

//CALCULATE MEAN DOY

var doycol = col_mosaic.select(['doy']).sort('system:time_start');

var mean_doy = doycol.map(function(img){
  var length = doycol.size();
  var list = doycol.toList(length);
  var index = list.indexOf(img);
  var previousIndex = ee.Algorithms.If(index.eq(0), index, index.subtract(1));
  var previousImage = ee.Image(list.get(previousIndex));
  var currentDate = ee.Date(previousImage.get('system:time_start'));
  var diffimg = doycol.filterDate(
                start, currentDate.advance(2, 'hour'));
  var reduce_image = diffimg.reduce(ee.Reducer.lastNonNull())
  .rename(['meanDOY']);
  // replace all masked values:
  var new_doy = img.add(reduce_image).divide(2).round().rename(['meanDOY']);
  return new_doy.addBands(reduce_image).rename(['meanDOY','DOYbefore']);
});

print('MEAN DOY', mean_doy);

// Incorporate Mean DOY bands into the collection 

var allcoll = col_mosaic.combine(mean_doy).sort('system:time_start');
var allcoll_fordiff = col_mosaic.sort('system:time_start');
print('all collection with mean doy bands', allcoll);

//CLOUD MASK EXCLUDING PIXELS WITH ST_B6 > 300
function maskL78sr(image) {
  // Bit 0 - Fill
  // Bit 1 - Dilated Cloud
  // Bit 2 - Unused
  // Bit 3 - Cloud
  // Bit 4 - Cloud Shadow
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // Apply the scaling factors to the appropriate bands.
  var opticalBands = image.select('Red','NIR','SWIR2').multiply(0.0000275).add(-0.2);
  var thermalBand = image.select('Thermal').multiply(0.00341802).add(149.0);
  var mascara = thermalBand.gte(300);
  var scal = image.addBands(opticalBands, null, true).addBands(thermalBand, null, true);
  var unm = scal.updateMask(qaMask).updateMask(saturationMask);
  var thermalmask = scal.mask(mascara);
  return unm.unmask(thermalmask);
}

// MAP THE CLOUD MASK FUNCTION OVER IMAGE COLLECTIONS

var maskedcollection_all = allcoll.map(maskL78sr).sort('system:time_start');
var maskedcollection_all_fordiff = allcoll_fordiff.map(maskL78sr).sort('system:time_start');

//test it (for L8)
var listOfImages_all = maskedcollection_all.toList(maskedcollection_all.size());
print('List of all images with cloud mask: ', listOfImages_all);
var cloud_masked_test = ee.Image(listOfImages_all.get(16));
Map.addLayer(cloud_masked_test, {bands: ['SWIR2', 'NIR', 'Red'], min: 0, max: 0.3}, 'Cloud masked - L8',false);

// CREATE NBR FUNCTION 
var nbr = function(img) {
 
    var index = img.normalizedDifference(['NIR', 'SWIR2'])
                  .select([0], ['NBR'])
                  .set('system:time_start', img.get('system:time_start'));
    return img.addBands(index).toFloat() ;
};

//Apply it
var ltCollection_all =  maskedcollection_all.map(nbr);
var ltCollection_all_fordiff =  maskedcollection_all_fordiff.map(nbr);

print('Filtered collection with NBR: ', ltCollection_all);

//Select only the bands of interest 
var select_nbr = ltCollection_all.select(['Red','NIR','SWIR2','Thermal','doy','NBR','meanDOY','DOYbefore']);

//Organize by ascendent date 
var select_nbr = select_nbr.sort('system:time_start');

print('Collection with selected bands: ', select_nbr);

//CALCULATE NBR DIFFERENCE: IMAGE(t) - IMAGE(t-1)
var nbrbands = ltCollection_all_fordiff.select(['NBR']).sort('system:time_start');

var NBRdiff = nbrbands.map(function(img){
  var length = nbrbands.size();
  var list = nbrbands.toList(length);
  var index = list.indexOf(img);
  var previousIndex = ee.Algorithms.If(index.eq(0), index, index.subtract(1));
  var previousImage = ee.Image(list.get(previousIndex));
  var currentDate = ee.Date(previousImage.get('system:time_start'));
  var diffimg = nbrbands.filterDate(
                start, currentDate.advance(2, 'hour'));
  var reduce_image = diffimg.reduce(ee.Reducer.lastNonNull())
  .rename(['NBRdiff']);
  // replace all masked values:
  var new_diff = img.subtract(reduce_image).rename(['NBRdiff']);
  return new_diff;
});

print('NBR diff',NBRdiff);

// Include NBRdiff collection into original collection
var allcoll_diff = select_nbr.combine(NBRdiff);
print('all collection with NBR and NBRdiff', allcoll_diff);

//PLOT TIME SERIES FOR SPECIFIC LAT-LONG
var point_18 = ee.Geometry.Point([-8.0323, 40.6588]);
Map.addLayer(point_18,{palette: ['FF8000']}, 'Point - time series chart',false);
var point_18 = ee.Feature(point_18);

// Reproject - the plot only works when using metric projections
var reproj = function(img) {
  var new_proj = img.reproject('EPSG:32629', null, 30);
  return new_proj;
};
var select_nbr_rep = select_nbr.map(reproj);

//NBR chart
var chart = ui.Chart.image.series({
    imageCollection: select_nbr_rep.select('NBR'),
    region: point_18.geometry()
    }).setOptions({
      interpolateNulls: true,
      lineWidth: 1,
      pointSize: 3,
      title: 'NBR time series before gap filling',
      vAxis: {title: 'Value'},
      hAxis: {title: 'Date', format: 'YYYY-MMM', gridlines: {count: 12}}
    });
print(chart);

var allcoll_diff_rep = allcoll_diff.map(reproj);

//NBR and deltaNBR chart
var chart = ui.Chart.image.series({
    imageCollection: allcoll_diff_rep.select('NBR','NBRdiff'),
    region: point_18.geometry()
    }).setOptions({
      interpolateNulls: true,
      lineWidth: 1,
      pointSize: 3,
      title: 'NBR and absolute NBR difference time series',
      vAxis: {title: 'Value'},
      hAxis: {title: 'Date', format: 'YYYY-MMM', gridlines: {count: 12}}
    });
print(chart);

//SWIR2, NIR and Red chart
var chart = ui.Chart.image.series({
    imageCollection: allcoll_diff_rep.select('SWIR2','NIR','Red'),
    region: point_18.geometry()
    }).setOptions({
      interpolateNulls: true,
      lineWidth: 1,
      pointSize: 3,
      title: 'RGB',
      vAxis: {title: 'Value'},
      hAxis: {title: 'Date', format: 'YYYY-MMM', gridlines: {count: 12}}
    });
print(chart);

//CREATE COMPOSITE IMAGE BASED ON MINIMUM deltaNBR (LARGEST NBR DROP)

var allcoll_diff_select = allcoll_diff.select(['NBRdiff','Red','NIR','SWIR2','Thermal','doy','NBR','meanDOY','DOYbefore']);

var min_NBRdiff = allcoll_diff_select.reduce(ee.Reducer.min(allcoll_diff_select.first().bandNames().size()))
  .rename(['NBRdiff','Red','NIR','SWIR2','Thermal','doy','NBR','meanDOY','DOYbefore']);

var min_NBRdiff = min_NBRdiff.select(['Red','NIR','SWIR2','Thermal','doy','NBR','NBRdiff','meanDOY','DOYbefore']);

print('Image reduced by min deltaNBR and other bands: ', min_NBRdiff);

var min_NBRdiff_pt = min_NBRdiff.clip(geom);

//Import and print burned areas (Annual Portuguese Fire Atlas)
var BA = ee.FeatureCollection('users/alanakneves/merge_2016_sup5');
Map.addLayer(BA, {color: 'FF0000'}, 'Burned area',false);

//Clip the minimum deltaNBR by the 2016 burned area limits
var min_NBRdiff_pt = min_NBRdiff_pt.clip(BA);

Map.addLayer(min_NBRdiff_pt, {bands: ['doy'], min: 0, max: 365}, 'DOY',false);
Map.addLayer(min_NBRdiff_pt, {bands: ['meanDOY'], min: 0, max: 365}, 'mean DOY');
Map.addLayer(min_NBRdiff_pt, {bands: ['DOYbefore'], min: 0, max: 365}, 'DOY before',false);

//Export the most relevant bands to Google Drive (delete /* and */ to run it)
/*
//NBRdiff
 Export.image.toDrive({
   image:min_NBRdiff_pt.select(['NBRdiff']).toFloat(),
   description: '2016_allPT_NBRdiff_yearcal_L7L8',
   folder: 'GEE-EXPORT',
   scale: 30,
   region: PT_int,
   maxPixels: 1e13,
   fileFormat: 'GeoTIFF',
   crs: 'EPSG:32629'
 });
//NBR
 Export.image.toDrive({
   image:min_NBRdiff_pt.select(['NBR']).toFloat(),
   description: '2016_allPT_NBR_yearcal_L7L8',
   folder: 'GEE-EXPORT',
   scale: 30,
   region: PT_int,
   maxPixels: 1e13,
   fileFormat: 'GeoTIFF',
   crs: 'EPSG:32629'
 }); 
//meanDOY
 Export.image.toDrive({
   image:min_NBRdiff_pt.select(['meanDOY']).toUint32(),
   description: '2016_allPT_meanDOY_yearcal_L7L8',
   folder: 'GEE-EXPORT',
   scale: 30,
   region: PT_int,
   maxPixels: 1e13,
   fileFormat: 'GeoTIFF',
   crs: 'EPSG:32629'
 });
 
//DOY
 Export.image.toDrive({
   image:min_NBRdiff_pt.select(['doy']).toUint32(),
   description: '2016_allPT_DOY_yearcal_L7L8',
   folder: 'GEE-EXPORT',
   scale: 30,
   region: PT_int,
   maxPixels: 1e13,
   fileFormat: 'GeoTIFF',
   crs: 'EPSG:32629'
 });
//DOYbefore
 Export.image.toDrive({
   image:min_NBRdiff_pt.select(['DOYbefore']).toUint32(),
   description: '2016_allPT_DOYbefore_yearcal_L7L8',
   folder: 'GEE-EXPORT',
   scale: 30,
   region: PT_int,
   maxPixels: 1e13,
   fileFormat: 'GeoTIFF',
   crs: 'EPSG:32629'
 });
*/
