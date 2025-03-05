var geom = ee.Geometry.Polygon(
        [[[-7.967291536150715, 41.40479767812003],
          [-7.967291536150715, 41.25939909190136],
          [-7.772284211931965, 41.25939909190136],
          [-7.772284211931965, 41.40479767812003]]], null, false);
Map.addLayer(geom, {color: 'yellow'}, 'ROI');

Map.centerObject(geom, 14);

// Alguns pontos de interesse (areas ardidas de 2020)
var campanho = ee.Geometry.MultiPoint([  [-7.9267835103483675, 41.323999310137985]]);
Map.addLayer(campanho, {color: 'purple'}, 'Campanho Point');

var outrospontosdeinteresse2020 = ee.Geometry.MultiPoint(
        [[-8.072464282674629, 41.54890664313104],
         [-8.592963657652168, 41.691533454745915],
         [-8.754657578341234, 41.83839985464384],
         [-7.113681399130716, 41.99828494740535]]);
Map.addLayer(outrospontosdeinteresse2020, {color: 'red'}, 'Outros pontos 2020');

// Areas ardidas ICNF 2020 do repositorio online (nao inclui a shp "outras areas ardidas" como fogo controlado, etc)
//Map.addLayer(ardida2020, {color: 'yellow', width: 2, fillColor: '00000000'},"Áreas Ardidas 2020");

///// 			INITIAL INPUT 				////

// YEAR TO BE PROCESSED
var currentyear = 2020;

// Data de inicio e fim para a recolha de imagens / analise
var start = ee.Date.fromYMD(currentyear, 1, 1).advance(-30, 'day'); // adicionado mes anterior à coleção
var end = ee.Date.fromYMD(currentyear, 12, 31);


// FILTER IMAGES BY PERIOD AND BBOX 
var filteredCollection_all_8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filter(ee.Filter.bounds(geom))
  .filterMetadata('CLOUD_COVER', 'less_than', 40)
  .filterDate(start, end).sort('system:time_start')
  // TEMPORARIO - APENAS PARA SELECIONAR PEQUENA SECCAO DAS IMAGENS PARA ACELERAR O VISUALIZADOR
  .map(function(image) { return image.clip(geom); }); 
   
var filteredCollection_all_7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filter(ee.Filter.bounds(geom))
  .filterMetadata('CLOUD_COVER', 'less_than', 40)  
  .filterDate(start, end).sort('system:time_start')
   // TEMPORARIO - APENAS PARA SELECIONAR PEQUENA SECCAO DAS IMAGENS PARA ACELERAR O VISUALIZADOR
  .map(function(image) { return image.clip(geom); }); 
  
//print('Filtered collection L8  - all year: ', filteredCollection_all_8);
//print('Filtered collection L7  - all year: ', filteredCollection_all_7);


///////////////////////  Kernel Gap Filling for L7 Collection : SLC-OFF  ///////////////////////
 //MC: outra abordagem: Stripe Error Correction for Landsat-7 Using Deep Learning: https://link.springer.com/article/10.1007/s41064-024-00306-x; https://github.com/ynsemrevrl/eliminating-stripe-errors

//input: filteredCollection_all_7
//output: l7_filtered_kernel

// Function to apply kernel gap filling to each image
function kernelGapFill(image) {
  // Define structuring elements
  var se3 = ee.Kernel.square(2);  // 3x3
  var se5 = ee.Kernel.square(4);  // 5x5
  var diag_kernel1 = ee.Kernel.fixed(3, 3, [[0, 0, 1], [0, 1, 0], [1, 0, 0]]);
  var diag_kernel2 = ee.Kernel.fixed(3, 3, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);

  // Apply Alternating Sequential Filter
  var eroded1 = image.focalMin({kernel: se3, iterations: 9}); 
  var dilated1 = eroded1.focalMax({kernel: se3, iterations: 9});

  var dilated2 = image.focalMax({kernel: se3, iterations: 9}); 
  var eroded2 = dilated2.focalMin({kernel: se3, iterations: 9});

  var asf_small = dilated1.add(eroded2).divide(2);

  // Apply diagonal median filtering
  var diag_filled1 = asf_small.focalMedian({kernel: diag_kernel1, iterations: 1});
  var diag_filled2 = asf_small.focalMedian({kernel: diag_kernel2, iterations: 1});

  var asf_diagonal = diag_filled1.add(diag_filled2).divide(2);

  // Final pass with larger structuring element
  var eroded_final = asf_diagonal.focalMin({kernel: se5, iterations: 9});
  var smoothed_image = eroded_final.focalMax({kernel: se5, iterations: 9});

  return smoothed_image.set('system:time_start', image.get('system:time_start'));
}

// Apply kernel gap filling to the collection
var l7_filtered = filteredCollection_all_7.map(kernelGapFill);
//print('colecaoL7_NBR_kernel', colecaoL7_NBR_kernel);

// Visualization of the first image after kernel gap filling
//Map.addLayer(l7_filtered.first(), {bands: ['NBR'], min: -1, max: 1}, 'colecaoL7_NBR_kernel');


  // --- APPLY SLC-OFF MASK AND BLEND INTERPOLATED DATA WITH ORIGINAL IMAGE ---
  
  // Function to apply SLC-off mask and blend with original image
  function applySLCOffMask(image) {
    var original_image = filteredCollection_all_7.filter(ee.Filter.eq('system:time_start', image.get('system:time_start'))).first();
  
    // Extract SLC-off mask from QA_PIXEL (Bit 0)
    var slc_off_mask = original_image.select('QA_PIXEL').toInt().bitwiseAnd(1).eq(1);
    var valid_pixel_mask = slc_off_mask.not();
  
    // Blend the valid pixels from the original image into the interpolated image //MC: de facto parece fazer mais sentido usar o where do que simplesmente o blend
    var final_image = image.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','ST_B6', 'SR_B7','QA_RADSAT', 'QA_PIXEL'])
                           .where(valid_pixel_mask, original_image.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','ST_B6', 'SR_B7', 'QA_RADSAT', 'QA_PIXEL']))
                           .set('system:time_start', image.get('system:time_start'));
  
    return final_image;
  }
  
  // Apply SLC-off mask and blending to the gap-filled collection
  var l7_filtered_kernel = l7_filtered.map(applySLCOffMask);
  //print('colecaoL7_NBR_kernelFilledCollection', l7_filtered_kernel);
  
  // Visualization of the first image after applying SLC-off mask and blending
  //Map.addLayer(l7_filtered_kernel.first(), {bands: ['SR_B7', 'SR_B4', 'SR_B3']}, 'colecaoL7_NBR_kernelFilledCollection');
  
  
/////////////////////// ADD DOY BAND AND RENAME L7 AND L8 BANDS ///////////////////////
  
  // ---  CREATE DATE BAND (DOY) FOR EACH IMAGE IN COLLECTION
  
  var addDate = function(image){
	  // Compute Day of Year (DOY)
    var doy = image.date().getRelative('day', 'year'); 
    var doyBand = ee.Image.constant(doy).uint16().add(1).rename('doy');
	// Compute Year
    var year = image.date().get('year');
    var yearBand = ee.Image.constant(year).uint16().rename('year');
     // Apply the original mask from SR_B4
    var mask = image.select('SR_B4').mask();
    doyBand = doyBand.updateMask(mask);
    yearBand = yearBand.updateMask(mask);
    // Add both
    return image.addBands([doyBand, yearBand]);
  };
  
  //apply it
  var withDate_l7 = l7_filtered_kernel.map(addDate).sort('system:time_start');
  //print('L7 Collection with doy band: ', withDate_l7);
  
  var withDate_l8 = filteredCollection_all_8.map(addDate).sort('system:time_start');
  //print('L8 Collection with doy band: ', withDate_l8);
  
  
  
  // --- RENAME L7 AND L5 BANDS TO MERGE THEM IN A UNIQUE COLLECTION
  
  //L7
  var rename_l7 = function(img) {
      var select = img.select('SR_B3','SR_B4','SR_B7','ST_B6','doy','year','QA_PIXEL','QA_RADSAT') //add band YEAR
                      .rename(['Red','NIR','SWIR2','Thermal','doy','year','QA_PIXEL','QA_RADSAT'])
                    .set('system:time_start', img.get('system:time_start'));
      return select;
  };
  
  var l7_renamed = withDate_l7.map(rename_l7);
  //print('L7 renamed:', l7_renamed);
  
  //L8
  var rename_l8 = function(img) {
      var select = img.select('SR_B4','SR_B5','SR_B7','ST_B10','doy','year','QA_PIXEL','QA_RADSAT')
                      .rename(['Red','NIR','SWIR2','Thermal','doy','year','QA_PIXEL','QA_RADSAT'])
                    .set('system:time_start', img.get('system:time_start'));
      return select;
  };
  
  var l8_renamed = withDate_l8.map(rename_l8);
  //print('L8 renamed:', l8_renamed);
  
  
///////////////////////  MERGE L7 AND L8  ///////////////////////
// output: col_mosaic 
  
  var filteredCollection_all = l8_renamed.merge(l7_renamed).sort('system:time_start');

  //print('Filtered collection - L7 and L8: ', filteredCollection_all);
  
   // ---  CREATE MOSAICS WITH IMAGES FROM THE SAME DATE
  
  function mosaicByDate(imcol){
    // imcol: An image collection
    // returns: An image collection
    var imlist = imcol.toList(imcol.size());
  
    var unique_dates = imlist.map(function(im){
      return ee.Image(im).date().format("YYYY-MM-dd"); // date: Returns the acquisition time of an image as a Date object. This helper function is equivalent to ee.Date(image.get('system:time_start')).
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
  
  
///////////////////////  MASCARA DE QUALIDADE baseada em QA_PIXEL e QA_RADSAT, e thermalBand > 300 ///////////////////////
// output: maskedcollection_all
  
  //CLOUD MASK EXCLUDING PIXELS WITH ST_B6 > 300
  function maskL78sr(image) {
    // Bit 0 - Fill
    // Bit 1 - Dilated Cloud
    // Bit 2 - Unused
    // Bit 3 - Cloud
    // Bit 4 - Cloud Shadow
  var qaMask = image.select('QA_PIXEL').toInt().bitwiseAnd(parseInt('11111', 2)).eq(0);
    var saturationMask = image.select('QA_RADSAT').eq(0);
  
    // Apply the scaling factors to the appropriate bands.
    // MC: fazer estas alterações no GEE vai converter inteiros em float e tornar o output mais pesado; uma solução é voltar a converter para inteiros 0-10000 por exemplo
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
  //print('maskedcollection_all',maskedcollection_all)
  //var maskedcollection_all_fordiff = allcoll_fordiff.map(maskL78sr).sort('system:time_start');
  
  
  ///////////////////////  ADD NBR BAND (NO TEMPORAL DIFFERENCES) ///////////////////////
// output: ltCollection_all

  // CREATE NBR FUNCTION 
  var nbr = function(img) {
   
      var index = img.normalizedDifference(['NIR', 'SWIR2'])
                    .select([0], ['NBR'])
                    .set('system:time_start', img.get('system:time_start'));
      return img.addBands(index).toFloat() ; // MC: aqui também é float
  };
  
  //Apply it
  var ltCollection_all =  maskedcollection_all.map(nbr); 
  
print('Filtered collection with NBR: ', ltCollection_all);


  
 
/////////////////////// FILTRO DE HAMPEL ///////////////////////

// input: ltCollection_all (image collection) com banda NBR
// OUTPUT: ltCollection_all_hampel: image collection, bandas ['Red','NIR','SWIR2','doy','NBR','year','Redb','NIRb','SWIR2b','doyb','NBRb','NBRdiff']

	// --- CALCULA NBRmedian  ---------
	// OUTPUT: IC NBRmedian

  //select only NBR bands
// a banda 'NBR' está em  ltCollection_all: porque é necessário calcular a IC NBRdiff
  var nbrbands = ltCollection_all.select(['NBR']);  // MC: não a dif de NBR , pois é a banda 'NBR' de NBRdiff, não e banda 'NBRdiff' // LL substitui NBRdiff por ltCollection_all

// MC: calcula a mediana de 3 sucessivas imagens?
  var NBRmedian = nbrbands.map(function(img){
    var length = nbrbands.size();
    var list = nbrbands.toList(length);
    var index = list.indexOf(img);
    //previous image
    var previousIndex = ee.Algorithms.If(index.eq(0), index, index.subtract(1));
    var previousImage = ee.Image(list.get(previousIndex));
    var currentDate = ee.Date(previousImage.get('system:time_start'));
    var diffimg = nbrbands.filterDate(
                  start, currentDate.advance(2, 'hour'));  // nbrbands é uma IC, por isso diffimg é também IC
    var reduce_image = diffimg.mosaic(); // o output de masaic é uma imagem
    //after image
    var nextIndex = ee.Algorithms.If(index.lt(length.subtract(1)), index.add(1), index); 
    var nextImage = ee.Image(list.get(nextIndex));
    var nextDate = ee.Date(nextImage.get('system:time_start'));
    var nextimg = nbrbands.filterDate(nextDate.advance(-2, 'hour'), end); //IC
    var nextimgrev = nextimg.sort('system:time_start', false);
    var reducenxt_image = nextimgrev.mosaic(); // devolve imagem
    // replace all masked values:
    var allthree = ee.ImageCollection([reduce_image, img, reducenxt_image]);
    var medianImage = allthree.reduce(ee.Reducer.median());
    medianImage = medianImage.updateMask(img.mask()).rename(['NBRmedian']);

  // Copia as propriedades originais (para preservar system:index e system:time_start)
  return medianImage.copyProperties(img, ['system:index', 'system:time_start']);
  });
  
  //print('NBRmedian',NBRmedian);

	// ---  CALCULA NBR-NBRmedian -----
	// output: IC diff_nbr_median

  //compute the difference between NBR and NBRmedian
  var diff_nbr_median = nbrbands.map(function(img){
    var length = nbrbands.size();
    var list = nbrbands.toList(length);
    var index = list.indexOf(img);
    var list_m = NBRmedian.toList(length);
    var image_m = ee.Image(list_m.get(index));
    var diffImage = img.subtract(image_m).rename(['NBRmedian_diff']);
  return diffImage.copyProperties(img, ['system:index', 'system:time_start']);
  });
  
  //print('diff_nbr_median',diff_nbr_median);


	// ---  NBR-NBRmedian (continuação) -----
	// output: imagens perc_25, etc (uma por ano)
	
  var perc_25 = diff_nbr_median.reduce(ee.Reducer.percentile([25]));
  var perc_75 = diff_nbr_median.reduce(ee.Reducer.percentile([75]));
  var iqr = perc_75.subtract(perc_25).multiply(1.5);
  var limiar_inf = perc_25.subtract(iqr);
  var limiar_sup = perc_75.add(iqr);
  
  //print('iqr',iqr);
  //Map.addLayer(iqr, null, 'iqr',false);

	// ---  substitui outliers em nbrbands por NBRmedium -----
// nota: nbrbands = NBRdiff.select(['NBR']), mas é o NBR, não a diferença
// output: replaced: IC com uma única banda

  var replace_iqr = function(img) {
    var length = nbrbands.size().add(10);
    var list = nbrbands.toList(length);
    var index = list.indexOf(img);
    
     // Garantir que o índice seja válido
    var validIndex = ee.Number(index).min(length.subtract(1)).max(0); // Índice entre 0 e (length - 1)
    var list_rep = NBRmedian.toList(length);
    var img_rep = ee.Image(list_rep.get(validIndex));
    
	// Substituição dos valores que excedem os limites
	var img_final = img.where(img.gt(limiar_sup), img_rep.unmask(img)).unmask(img);
	img_final = img_final.where(img.lt(limiar_inf), img_rep.unmask(img)).unmask(img);
	img_final = img_final.rename(['NBRrep']);

  return img_final.copyProperties(img, ['system:index', 'system:time_start']);  
  };
  
  
  var replaced = nbrbands.map(replace_iqr); // NBRrep são os valores de nbrbands, e por isso de NBR, e não de diferença de NBR
  print('replaced',replaced);
  
  // Combina as bandas processadas com a coleção original
// Nesta abordagem, para cada imagem original é adicionado o NBRmedian, NBRmedian_diff e NBRrep
var ltCollection_all_hampel = ltCollection_all.map(function(img) {
  var id = img.get('system:index');  // Preserva o ID original

  // Busca as imagens correspondentes nas coleções processadas com base em system:index
  var medianImg = ee.Image(NBRmedian.filter(ee.Filter.eq('system:index', id)).first());
  var diffImg = ee.Image(diff_nbr_median.filter(ee.Filter.eq('system:index', id)).first());
  var repImg = ee.Image(replaced.filter(ee.Filter.eq('system:index', id)).first());

  // Caso alguma imagem não seja encontrada, substitui por uma imagem constante (0)
  medianImg = ee.Algorithms.If(medianImg, medianImg, ee.Image.constant(0).rename('NBRmedian'));
  diffImg = ee.Algorithms.If(diffImg, diffImg, ee.Image.constant(0).rename('NBRmedian_diff'));
  repImg = ee.Algorithms.If(repImg, repImg, ee.Image.constant(0).rename('NBRrep'));

  return img.addBands(ee.Image(medianImg))
            .addBands(ee.Image(diffImg))
            .addBands(ee.Image(repImg));
}).sort('system:time_start');

print('Fixed Hampel Collection:', ltCollection_all_hampel);
  

	///////////  R_NBR ////////////

	//CALCULAR o ÍNDICE R_NBR que depois é usado para o REDUCER  
	


 //Select only the bands of interest  
  var select_bands_hampel = ltCollection_all_hampel.select(['Red', 'NIR', 'SWIR2', 'doy', 'year', 'NBRrep']) // Colocado NBRrep ao inves de NBR
   
  //Organize by ascendent date 
  var sorted = select_bands_hampel.sort('system:time_start');
  
  // Convert to list
  var list = sorted.toList(sorted.size());

	// nova coleção iterando por índices
var withBefore = ee.ImageCollection(
  ee.List.sequence(0, list.size().subtract(1)).map(function(i) {
    i = ee.Number(i);
    
    // Imagem atual (forçada a float)
    var current = ee.Image(list.get(i)).toFloat();

    // Índice da imagem anterior
    var prevIndex = i.subtract(1);
    
    // Se existe anterior, usa-a e converte a float; caso contrário, dummy zeros float
    var prevImage = ee.Image(
      ee.Algorithms.If(
        prevIndex.gte(0),
        ee.Image(list.get(prevIndex))
          .select(['Red','NIR','SWIR2','doy','year','NBRrep'])
          .rename(['Redb','NIRb','SWIR2b','doyb','yearb','NBRrep_b'])
          .toFloat(),
        ee.Image.constant([0, 0, 0, 0, 0, 0])
          .rename(['Redb','NIRb','SWIR2b','doyb','yearb','NBRrep_b'])
          .toFloat()
      )
    );

    // Combina a atual com a anterior
    var combined = current.addBands(prevImage);

    // Copia as propriedades originais da imagem atual (por ex. time_start, index)
    // para manter metadados (opcional, mas recomendável)
    return combined.copyProperties(current, current.propertyNames());
  })
);

// Visualiza
//Map.addLayer(withBefore, {bands: ['SWIR2b','NIRb','Redb']}, 'NBRdiff_hampel');
print('NBRdiff_hampel', withBefore);



		// ---  Function to compute NBR change ratio  ------

	// Compute the ratio: (dNBR)/(NBR + 1.001)
	//(NBR_Prefire - NBR_Postfire) / (NBR_Postfire + 1.001)

var addNBRChange = function(image) {
	  // Ensure NBR comes only from the current year
    var currentNBR = image.updateMask(image.select('year').eq(currentyear));
      // Compute dNBR as the difference between the two NBR values
	var dNBR = image.select('NBRrep_b').subtract(image.select('NBRrep'));
    // Compute the ratio: (dNBR) / (NBR + 1.001) - The small constant (1.001) helps avoid division by zero
  var NBRChange = dNBR.divide(image.select('NBRrep').add(1.001))
                        .rename('NBRChange');
    // Add the new band to the image
  return image.addBands(NBRChange);
};

// Map the function over the NBRdiff collection to get the final result
var NBRdiff_NBRchange = withBefore.map(addNBRChange);


///////////////////// REDUCER PARA OBTER COMPOSITO ANUAL //////////////////////////

//REDUCE BY MAX_NBRCHANGE

// LL Substituí ee.reduce por qualityMosaic as "The qualityMosaic() function selects the image (per-pixel) with the HIGHEST quality-band-score to contribute to the resulting mosaic". That is, the values are taken from the image where NBRChange is maximum.

var MAX_NBRCHANGE = NBRdiff_NBRchange.qualityMosaic('NBRChange');

// Print and check the result
print('MAX_NBRCHANGE composite:', MAX_NBRCHANGE);


/////  VALIDATION OF THE LIMITATION   
// CALCULATE HOW MANY PIXELS WITHIN A CERTAIN AREA ARE FROM THE PREVIOUS YEAR
// SHOULD EXIST SOME BUT NOT MANY

// Create a binary mask where year == targetYear
var yearMaskcurrent = MAX_NBRCHANGE.clip(geom)
                            .select('year')
                            .eq(currentyear); // add -1 para o previous year

// Count non-masked pixels using reduceRegion
var pixelCountcurrent = yearMaskcurrent.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geom,  // Use the defined region of interest
    scale: 30,       // Landsat resolution
    maxPixels: 1e13
});

// Print the total number of pixels for the selected year
print('Total number of pixels for year', currentyear, pixelCountcurrent.get('year'));


// previous year
// Create a binary mask where year =! targetYear
var yearMaskprevious = MAX_NBRCHANGE.clip(geom)
                            .select('year')
                            .eq(currentyear-1); // add -1 para o previous year

// Count non-masked pixels using reduceRegion
var pixelCountprevious = yearMaskprevious.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geom,  // Use the defined region of interest
    scale: 30,       // Landsat resolution
    maxPixels: 1e13
});
print('Total number of pixels for previous year', currentyear - 1, pixelCountprevious.get('year'));

//SEMHAMPEL - Total number of pixels for year 2019 [GEOM area] - 12049.835294117634
//SEMHAMPEL - Total number of pixels for year 2020 [GEOM area] - 378423.8666666629
//SEMHAMPEL - Corresponde a cerca de 3%


///////////////////// HOMOGENEITY //////////////////////////

// Apply Haralick homogeneity based on DOY_Current
var doyCurrent = MAX_NBRCHANGE.select('doy').toUint16(); // Convert to integer as Image.glcmTexture: Only 32-bit or smaller integer types are currently supported.

var doyCurrentTexture = doyCurrent.glcmTexture({size: 3}).select('doy_idm').rename('doy_homogeneity'); 
  // das várias bandas criadas a de interesse é : IDM: f5, Inverse Difference Moment; measures the homogeneity

// Add the selected band to MAX_NBRCHANGE (currentComposite)
var AnualComposite = MAX_NBRCHANGE.addBands(doyCurrentTexture);
print('Annual Composite ' + currentyear ,AnualComposite);


//////      OUTPUT      //////

// OUTPUT: 'AnualComposite' com multibandas - 

// Clip para visualizar pequena area
var AnualComposite_geom = AnualComposite.clip(geom);
Map.addLayer(AnualComposite_geom, {bands: ['NBRChange']}, 'Annual Composite ' + currentyear + 'NBRChange');
Map.addLayer(AnualComposite_geom, {bands: ['SWIR2', 'NIR', 'Red']}, 'Annual Composite ' + currentyear + 'RGB');
Map.addLayer(AnualComposite_geom, {bands: ['doy_homogeneity']}, 'Annual Composite ' + currentyear + 'DOY Homogeneity');


//////      EXPORT      //////

//Export MULTIBAND RASTER
 Export.image.toDrive({
   image:AnualComposite_geom.toFloat(),
   description: 'Composito'+currentyear,
   folder: 'GEE_EXPORT_ATLAS',
   scale: 30,
   region: geom, //ajustar
   maxPixels: 1e13,
   fileFormat: 'GeoTIFF',
   crs: 'EPSG:4326' // -> crs WGS84 - ADAPTAR PARA EPSG:3763 (ETRS89 / Portugal TM06)?
 });