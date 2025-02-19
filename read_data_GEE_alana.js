//Visualization of time series
//replace the original NBR series by the NBR median in outliers
//uses the difference between NBR and NBRmedian

//Import Portugal limit
//Feature collection
var PT_int = ee.FeatureCollection('users/alanakneves/NSR_buffer');
//Collect the feature collection geometry
var geom = PT_int.geometry();

//print(geom);

Map.centerObject(geom, 6);

var start = ee.Date('2019-05-01');
var end = ee.Date('2019-12-31');

// FILTER IMAGES BY PERIOD AND BBOX 

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

// APPLY FOCAL MEAN TO FILL L7 GAPS

var filteredCollection_all_7 = filteredCollection_all_7.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','ST_B6','QA_PIXEL','QA_RADSAT']);

var filtered_focal_mean = function(img) {
  var img_fill_1 = img.select('QA_PIXEL','QA_RADSAT').focal_mode(1, 'square', 'pixels', 8);
  var img_fill_2 = img.select('SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','ST_B6').focal_mean(1, 'square', 'pixels', 8);
  var img_fill = img_fill_2.addBands(img_fill_1);
  var img_final = img_fill.blend(img);
  return img_final.int().copyProperties(img,['system:time_start']);
};

var l7_filtered = filteredCollection_all_7.map(filtered_focal_mean);

//print('L7 after focal mean', l7_filtered);

//CREATE DATE BAND (DOY) FOR EACH IMAGE IN COLLECTION

var addDate = function(image){
  var doy = image.date().getRelative('day', 'year');
  var doyBand = ee.Image.constant(doy).uint16().add(1).rename('doy');
  doyBand = doyBand.updateMask(image.select('SR_B4').mask());
  
  return image.addBands(doyBand);
};

//apply it
var withDate_l7 = l7_filtered.map(addDate).sort('system:time_start');
//print('L7 Collection with doy band: ', withDate_l7);

var withDate_l8 = filteredCollection_all_8.map(addDate).sort('system:time_start');
//print('L8 Collection with doy band: ', withDate_l8);

//RENAME L7 AND L5 BANDS TO MERGE THEM IN A UNIQUE COLLECTION

//L7
var rename_l7 = function(img) {
    var select = img.select('SR_B3','SR_B4','SR_B7','ST_B6','doy','QA_PIXEL','QA_RADSAT')
                    .rename(['Red','NIR','SWIR2','Thermal','doy','QA_PIXEL','QA_RADSAT'])
                  .set('system:time_start', img.get('system:time_start'));
    return select;
};

var l7_renamed = withDate_l7.map(rename_l7);
//print('L7 renamed:', l7_renamed);

//L8
var rename_l8 = function(img) {
    var select = img.select('SR_B4','SR_B5','SR_B7','ST_B10','doy','QA_PIXEL','QA_RADSAT')
                    .rename(['Red','NIR','SWIR2','Thermal','doy','QA_PIXEL','QA_RADSAT'])
                  .set('system:time_start', img.get('system:time_start'));
    return select;
};

var l8_renamed = withDate_l8.map(rename_l8);
//print('L8 renamed:', l8_renamed);

// MERGE L7 AND L8

var filteredCollection_all = l8_renamed.merge(l7_renamed).sort('system:time_start');
//print('Filtered collection - L7 and L8: ', filteredCollection_all);

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

//print('Mosaics of the same date', col_mosaic);

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

var maskedcollection_all = col_mosaic.map(maskL78sr).sort('system:time_start');
//var maskedcollection_all_fordiff = allcoll_fordiff.map(maskL78sr).sort('system:time_start');

// CREATE NBR FUNCTION 
var nbr = function(img) {
 
    var index = img.normalizedDifference(['NIR', 'SWIR2'])
                  .select([0], ['NBR'])
                  .set('system:time_start', img.get('system:time_start'));
    return img.addBands(index).toFloat() ;
};

//Apply it
var ltCollection_all =  maskedcollection_all.map(nbr);

//print('Filtered collection with NBR: ', ltCollection_all);

//Select only the bands of interest 
var select_bands = ltCollection_all.select(['Red','NIR','SWIR2','doy','NBR']);

//Organize by ascendent date 
var select_bands = select_bands.sort('system:time_start');

var NBRdiff = select_bands.map(function(img){
  var length = select_bands.size();
  var list = select_bands.toList(length);
  var index = list.indexOf(img);
  var previousIndex = ee.Algorithms.If(index.eq(0), index, index.subtract(1));
  var previousImage = ee.Image(list.get(previousIndex));
  var currentDate = ee.Date(previousImage.get('system:time_start'));
  var diffimg = select_bands.filterDate(
                start, currentDate.advance(2, 'hour'));
  var reduce_image = diffimg.reduce(ee.Reducer.lastNonNull())
  .rename(['Redb','NIRb','SWIR2b','doyb','NBRb']);
  // replace all masked values:
  var new_diff = img.select(['NBR']).subtract(reduce_image.select(['NBRb'])).rename(['NBRdiff']);
  var result = img.addBands(reduce_image).rename(['Red','NIR','SWIR2','doy','NBR','Redb','NIRb','SWIR2b','doyb','NBRb']);
  return result.addBands(new_diff).rename(['Red','NIR','SWIR2','doy','NBR','Redb','NIRb','SWIR2b','doyb','NBRb','NBRdiff']);
});

//print('NBR diff - all bands',NBRdiff);

//select only NBR bands
var nbrbands = NBRdiff.select(['NBR']);

var NBRmedian = nbrbands.map(function(img){
  var length = nbrbands.size();
  var list = nbrbands.toList(length);
  var index = list.indexOf(img);
  //previous image
  var previousIndex = ee.Algorithms.If(index.eq(0), index, index.subtract(1));
  var previousImage = ee.Image(list.get(previousIndex));
  var currentDate = ee.Date(previousImage.get('system:time_start'));
  var diffimg = nbrbands.filterDate(
                start, currentDate.advance(2, 'hour'));
  var reduce_image = diffimg.mosaic();
  //after image
  var nextIndex = ee.Algorithms.If(index.eq(81), index, index.add(1));
  var nextImage = ee.Image(list.get(nextIndex));
  var nextDate = ee.Date(nextImage.get('system:time_start'));
  var nextimg = nbrbands.filterDate(nextDate.advance(-2, 'hour'), end);
  var nextimgrev = nextimg.sort('system:time_start', false);
  var reducenxt_image = nextimgrev.mosaic();
  // replace all masked values:
  var allthree = ee.ImageCollection([reduce_image, img, reducenxt_image]);
  var medianImage = allthree.reduce(ee.Reducer.median());
  return medianImage.updateMask(img.mask()).rename(['NBRmedian']);
});

//print('NBRmedian',NBRmedian);

//compute the difference between NBR and NBRmedian
var diff_nbr_median = nbrbands.map(function(img){
  var length = nbrbands.size();
  var list = nbrbands.toList(length);
  var index = list.indexOf(img);
  var list_m = NBRmedian.toList(length);
  var Image_m = ee.Image(list_m.get(index));
  return img.subtract(Image_m).rename(['NBRmedian_diff']);
});

//print('diff_nbr_median',diff_nbr_median);

var perc_25 = diff_nbr_median.reduce(ee.Reducer.percentile([25]));
var perc_75 = diff_nbr_median.reduce(ee.Reducer.percentile([75]));
var iqr = perc_75.subtract(perc_25);
var iqr = iqr.multiply(1.5);
var limiar_inf = perc_25.subtract(iqr);
var limiar_sup = perc_75.add(iqr);

//print('iqr',iqr);
//Map.addLayer(iqr, null, 'iqr',false);

var replace_iqr = function(img) {
  var length = nbrbands.size();
  var list = nbrbands.toList(length);
  var index = list.indexOf(img);
  var list_rep = NBRmedian.toList(length);
  var img_rep = ee.Image(list_rep.get(index));
  var img_final = img.where(img.gt(limiar_sup), img_rep);
  var img_final = img_final.where(img.lt(limiar_inf), img_rep);
  return img_final.rename(['NBRrep']);
};

var replaced = nbrbands.map(replace_iqr);
//print('replaced',replaced);

var allcoll_rep = NBRdiff.combine(NBRmedian);
var allcoll_rep = allcoll_rep.combine(diff_nbr_median);
var allcoll_rep = allcoll_rep.combine(replaced).sort('system:time_start');
print('replaced',allcoll_rep);


/*
var list_allcoll_diff = allcoll_diff.sort('system:time_start').toList(allcoll_diff.size());
//print('List of all NBRdiff images', list_allcoll_diff);
var img_test = ee.Image(list_allcoll_diff.get(32));
Map.addLayer(img_test, {bands: ['SWIR2', 'NIR', 'Red'], min: 0, max: 0.3}, '14/jan',false);
var img_test = ee.Image(list_allcoll_diff.get(9));
Map.addLayer(img_test, {bands: ['SWIR2', 'NIR', 'Red'], min: 0, max: 0.3}, '22/01',false);
var img_test = ee.Image(list_allcoll_diff.get(14));
Map.addLayer(img_test, {bands: ['SWIR2', 'NIR', 'Red'], min: 0, max: 0.3}, '7/fev',false);
*/
/*
var point_18 = ee.Geometry.Point([37.4679, -12.4057]);
Map.addLayer(point_18,{palette: ['FF8000']}, 'Ponto');
var point_18 = ee.Feature(point_18);

// the plot only works when using metric projections
var reproj = function(img) {
  var new_proj = img.reproject('EPSG:32629', null, 30);
  return new_proj;
};
var select_nbr_rep = allcoll_rep.map(reproj);

var chart = ui.Chart.image.series({
    imageCollection: select_nbr_rep.select('NBRmedian','NBRrep','NBRmedian_diff','NBR'),
    region: point_18.geometry()
    }).setOptions({
      interpolateNulls: true,
      lineWidth: 1,
      pointSize: 3,
      title: 'NBR int time series',
      vAxis: {title: 'Value'},
      hAxis: {title: 'Date', format: 'YYYY-MMM', gridlines: {count: 12}}
    });
print(chart);
*/

//REDUCE BY MIN NBR
var allcoll_diff_select = allcoll_rep.select(['NBRrep','Red','NIR','SWIR2','doy','Redb','NIRb','SWIR2b','doyb']);

var min_NBRdiff = allcoll_diff_select.reduce(ee.Reducer.min(allcoll_diff_select.first().bandNames().size()))
  .rename(['NBRrep','Red','NIR','SWIR2','doy','Redb','NIRb','SWIR2b','doyb']);

print('Image reduced by min NBR and other bands: ', min_NBRdiff);

var min_NBRdiff_pt = min_NBRdiff.clip(geom);

Map.addLayer(min_NBRdiff_pt, {bands: ['SWIR2', 'NIR', 'Red'], min: 0, max: 0.3}, 'min NBR',false);

/*
var studyRegion = ee.Geometry.Rectangle([
   [37.105, -12.5826],
   [37.4731, -12.8451]
]);

var studyRegion = ee.FeatureCollection(ee.Feature(studyRegion));
*/

//var studyRegion = ee.FeatureCollection('users/alanakneves/NSR_smallarea_BA');

var styling = {color: 'yellow', fillColor: '00000000'};
Map.addLayer(PT_int.style(styling),null,'area');
//Map.centerObject(studyRegion, 10.0);

/*
//Export blended image
 Export.image.toDrive({
   image:min_NBRdiff_pt.toFloat(),
   description: 'NSR_2019_L7L8_minNBR_wbuffer_40cc',
   folder: 'GEE-EXPORT_DL',
   scale: 30,
   region: PT_int,
   maxPixels: 1e13,
   fileFormat: 'GeoTIFF',
   crs: 'EPSG:4326'
 });
*/
