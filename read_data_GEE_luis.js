var geom = 
    /* color: #98ff00 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-7.967291536150715, 41.40479767812003],
          [-7.967291536150715, 41.25939909190136],
          [-7.772284211931965, 41.25939909190136],
          [-7.772284211931965, 41.40479767812003]]], null, false),
    campanho = /* color: #d63000 */ee.Geometry.MultiPoint(
        [[-7.9267835103483675, 41.323999310137985],  /// ponto output
         [-7.154017760095339, 41.41880931678191]]),
    cab_bastos_Set2020 = /* color: #d63000 */ee.Geometry.MultiPoint(
        [[-8.072464282674629, 41.54890664313104],
         [-8.592963657652168, 41.691533454745915],
         [-8.754657578341234, 41.83839985464384],
         [-7.113681399130716, 41.99828494740535]]);

// Areas ardidas ICNF 2020 do repositorio online (nao inclui a shp "outras areas ardidas" como fogo controlado, etc)
//Map.addLayer(ardida2020, {color: 'yellow', width: 2, fillColor: '00000000'},"Áreas Ardidas 2020");

Map.centerObject(geom, 11);


// Data de inicio e fim para a recolha de imagens / analise
var start = ee.Date('2019-01-01');
var end = ee.Date('2020-12-31');

// FILTER IMAGES BY PERIOD AND BBOX 

var filteredCollection_all_8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filter(ee.Filter.bounds(geom))
  .filterMetadata('CLOUD_COVER', 'less_than', 40)
  .filterDate(start, end).sort('system:time_start');
   
var filteredCollection_all_7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filter(ee.Filter.bounds(geom))
  .filterMetadata('CLOUD_COVER', 'less_than', 40)  // MC: demasiado restrito?
  .filterDate(start, end).sort('system:time_start');
  
print('Filtered collection L8  - all year: ', filteredCollection_all_8);
print('Filtered collection L7  - all year: ', filteredCollection_all_7);


/* 

MC: outra abordagem: Stripe Error Correction for Landsat-7 Using Deep Learning: https://link.springer.com/article/10.1007/s41064-024-00306-x; https://github.com/ynsemrevrl/eliminating-stripe-errors

// FOCAL MEAN ERA A ABORDAGEM DA ALANA PARA PREENCHER AS FALHAS DO L7 EM TERMOS DE SLC-OFF 
// ESTA ABORDAGEM FOI SUBSTITUIDA POR UM FILTRO KERNEL

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

*/


///////////////////////  Kernel Gap Filling for Entire Collection : SLC-OFF  ///////////////////////

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
  print('colecaoL7_NBR_kernelFilledCollection', l7_filtered_kernel);
  
  // Visualization of the first image after applying SLC-off mask and blending
  Map.addLayer(l7_filtered_kernel.first(), {bands: ['SR_B7', 'SR_B4', 'SR_B3']}, 'colecaoL7_NBR_kernelFilledCollection');
  
/////////////////////// ADD DOY BAND AND RENAME L7 AND L8 BANDS ///////////////////////
  
  //CREATE DATE BAND (DOY) FOR EACH IMAGE IN COLLECTION
  
  var addDate = function(image){
    var doy = image.date().getRelative('day', 'year');
    var doyBand = ee.Image.constant(doy).uint16().add(1).rename('doy');
    doyBand = doyBand.updateMask(image.select('SR_B4').mask());
    
    return image.addBands(doyBand);
  };
  
  //apply it
  var withDate_l7 = l7_filtered_kernel.map(addDate).sort('system:time_start');
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
  
  
///////////////////////  MERGE L7 AND L8  ///////////////////////
// output: col_mosaic 
  
  var filteredCollection_all = l8_renamed.merge(l7_renamed).sort('system:time_start');

  //print('Filtered collection - L7 and L8: ', filteredCollection_all);
  
  // CREATE MOSAICS WITH IMAGES FROM THE SAME DATE
  
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
    var mascara = thermalBand.gte(300);  // MC: seria melhor ter estas constantes definidas no cabeçalho
    var scal = image.addBands(opticalBands, null, true).addBands(thermalBand, null, true);
    var unm = scal.updateMask(qaMask).updateMask(saturationMask);
    var thermalmask = scal.mask(mascara);
    return unm.unmask(thermalmask);
  }
  
  // MAP THE CLOUD MASK FUNCTION OVER IMAGE COLLECTIONS
  
  var maskedcollection_all = col_mosaic.map(maskL78sr).sort('system:time_start');
  print('maskedcollection_all',maskedcollection_all)
  //var maskedcollection_all_fordiff = allcoll_fordiff.map(maskL78sr).sort('system:time_start');
  
  
  ///////////////////////  ADD NBR BAND (NO TEMPORAL DIFFERENCES) ///////////////////////
// output: ltCollection_all
// output: ltCollection_all_prefiltro // image collection com bandas ['NBR', 'Red','NIR','SWIR2','doy']

  // CREATE NBR FUNCTION 
  var nbr = function(img) {
   
      var index = img.normalizedDifference(['NIR', 'SWIR2'])
                    .select([0], ['NBR'])
                    .set('system:time_start', img.get('system:time_start'));
      return img.addBands(index).toFloat() ; // MC: aqui também é float
  };
  
  //Apply it
  var ltCollection_all =  maskedcollection_all.map(nbr); // só tem 5 bandas: nbr, red, nir, swir2 e doy

print('Filtered collection with NBR: ', ltCollection_all);

// MC: Novo
  // .REDUCE BY MIN NBR AQUI APENAS PARA SE TER IMAGEM PREFILTRO HAMPEL
  // PARA SER MAIS LEVE, ESTA APENAS PARA A AREA "geom"
  
  var ltCollection_all_prefiltro = ltCollection_all.select(['NBR', 'Red','NIR','SWIR2','doy']);
  print('Filtered collection with NBR (prefiltro) ', ltCollection_all_prefiltro);

/////////////////////// REDUCE COLLECTION TO IMAGE WITH MIN(NBR) /////////////////////// ltCollection_all_prefiltro_min_NBRdiff: Não é usado fora deste bloco
// OUTPUT: ltCollection_all_prefiltro_min_NBRdiff_GEOM 
// MC: NBRdiff é enganador

// aplicar critério min NBR para redução
  var ltCollection_all_prefiltro_min_NBRdiff = ltCollection_all_prefiltro.reduce(ee.Reducer.min(ltCollection_all_prefiltro.first().bandNames().size()))
    .rename(['NBR', 'Red','NIR','SWIR2','doy']);
  
  print('Image reduced by min NBR and other bands: ', ltCollection_all_prefiltro_min_NBRdiff);
  
  var ltCollection_all_prefiltro_min_NBRdiff_GEOM = ltCollection_all_prefiltro_min_NBRdiff.clip(geom);
  
  Map.addLayer(ltCollection_all_prefiltro_min_NBRdiff_GEOM, {bands: ['SWIR2', 'NIR', 'Red']}, 'min NBR_RGB_PREFILTRO',false);
  
  Map.addLayer(_GEOM, {bands: ['NBR']}, 'min NBR prefiltro',false);
  
////////////////////////////////////////////////////////////////  
/////////////////////// FILTRO DE HAMPEL ///////////////////////
////////////////////////////////////////////////////////////////
// outputs: IC NBRdiff; IC NBRmedian; IC diff_nbr_median ; IC allcol_rep 
// É mesmo NBRdiff que se quer?

/////////////////////// A. CALCULA NBRdiff COMO ALANA ///////////////////////
// input: ltCollection_all (image collection) com banda NBR
// OUTPUT: NBRdiff: image collection, bandas ['Red','NIR','SWIR2','doy','NBR','Redb','NIRb','SWIR2b','doyb','NBRb','NBRdiff']

  //Select only the bands of interest 
  var select_bands = ltCollection_all.select(['Red','NIR','SWIR2','doy','NBR']); // MC: são as unicas bandas desta IC
  
  //Organize by ascendent date 
  var select_bands = select_bands.sort('system:time_start');

// a 1a banda de NBRdiff chamada 'NBRdiff' vai ser de facto a diferença dos NBR
  var NBRdiff = select_bands.map(function(img){
    var length = select_bands.size();
    var list = select_bands.toList(length);
    var index = list.indexOf(img);
    var previousIndex = ee.Algorithms.If(index.eq(0), index, index.subtract(1));     
    var previousImage = ee.Image(list.get(previousIndex));
    var currentDate = ee.Date(previousImage.get('system:time_start'));
    var diffimg = select_bands.filterDate(
                  start, currentDate.advance(2, 'hour')); //IC
    var reduce_image = diffimg.reduce(ee.Reducer.lastNonNull()) // image
    .rename(['Redb','NIRb','SWIR2b','doyb','NBRb']);
    // replace all masked values:
    var new_diff = img.select(['NBR']).subtract(reduce_image.select(['NBRb'])).rename(['NBRdiff']);
    var result = img.addBands(reduce_image).rename(['Red','NIR','SWIR2','doy','NBR','Redb','NIRb','SWIR2b','doyb','NBRb']);
    return result.addBands(new_diff).rename(['Red','NIR','SWIR2','doy','NBR','Redb','NIRb','SWIR2b','doyb','NBRb','NBRdiff']);
  });
  
  print('NBR diff - all bands',NBRdiff);

///////////////////// CALCULA NBRmedian COMO ALANA //////////////////////////
// OUTPUT: IC NBRmedian

  //select only NBR bands
// a banda 'NBR' está em  ltCollection_all: porque é necessário calcular a IC NBRdiff
  var nbrbands = NBRdiff.select(['NBR']);  // MC: não a dif de NBR , pois é a banda 'NBR' de NBRdiff, não e banda 'NBRdiff'

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
        //  NESTE SEQ. DE CODIGO DA ALANA É NECESSARIO AJUSTAR O NUMERO DE index.eq(X) CONSOANTE O
        //   NUMERO DE IMAGENS EXISTENTES NA COLECAO 
    // MC: para resolver isso, fazer   var nextIndex = ee.Algorithms.If(index.lt(length), index, index.add(1)); //
    var nextIndex = ee.Algorithms.If(index.eq(66), index, index.add(1)); // MC: colocar length em vez de 66
    var nextImage = ee.Image(list.get(nextIndex));
    var nextDate = ee.Date(nextImage.get('system:time_start'));
    var nextimg = nbrbands.filterDate(nextDate.advance(-2, 'hour'), end); //IC
    var nextimgrev = nextimg.sort('system:time_start', false);
    var reducenxt_image = nextimgrev.mosaic(); // devolve imagem
    // replace all masked values:
    var allthree = ee.ImageCollection([reduce_image, img, reducenxt_image]);
    var medianImage = allthree.reduce(ee.Reducer.median());
    return medianImage.updateMask(img.mask()).rename(['NBRmedian']);
  });
  
  //print('NBRmedian',NBRmedian);

///////////////////// CALCULA NBR-NBRmedian COMO ALANA //////////////////////////
// output: IC diff_nbr_median

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


///////////////////// CALCULA NBR-NBRmedian COMO ALANA //////////////////////////
// output: imagens perc_25, etc (uma por ano)
  var perc_25 = diff_nbr_median.reduce(ee.Reducer.percentile([25]));
  var perc_75 = diff_nbr_median.reduce(ee.Reducer.percentile([75]));
  var iqr = perc_75.subtract(perc_25);
  var iqr = iqr.multiply(1.5);
  var limiar_inf = perc_25.subtract(iqr);
  var limiar_sup = perc_75.add(iqr);
  
  //print('iqr',iqr);
  //Map.addLayer(iqr, null, 'iqr',false);

///////////////////// substitui outliers em nbrbands por NBRmedium COMO ALANA //////////////////////////
// nota: nbrbands = NBRdiff.select(['NBR']), mas é o NBR, não a diferença
// output: replaced: IC com uma única banda

  var replace_iqr = function(img) {
    var length = nbrbands.size().add(10);
    var list = nbrbands.toList(length);
    var index = list.indexOf(img);
    
     // Garantir que o índice seja válido
    var validIndex = ee.Number(index).min(length.add(2)).max(0); // Índice entre 0 e (length - 1)
    var list_rep = NBRmedian.toList(length);
    
    var img_rep = ee.Image(list_rep.get(index));
    var img_final = img.where(img.gt(limiar_sup), img_rep);
    var img_final = img_final.where(img.lt(limiar_inf), img_rep);
    return img_final.rename(['NBRrep']);  
  };
  
  var replaced = nbrbands.map(replace_iqr); // NBRrep são os valores de nbrbands, e por isso de NBR, e não de diferença de NBR
  print('replaced',replaced);
  
  var allcoll_rep = NBRdiff.combine(NBRmedian); // combines Makes a new collection that is a copy of the images in primary, adding all the bands from the image in secondary with a matching ID.
  var allcoll_rep = allcoll_rep.combine(diff_nbr_median);
  var allcoll_rep = allcoll_rep.combine(replaced).sort('system:time_start');
  print('replaced',allcoll_rep);
  
///////////////////// FINALMENTE, FAZ REDUCE PARA OBTER OUTPUT PARA EXPORTAR //////////////////////////
// o reducer é o mínimo de NBRrep, que são os valores de NBR corrigidos para os outliers
// parece que não usa a banda 'NBRdiff'
//REDUCE BY MIN NBR
var allcoll_diff_select = allcoll_rep.select(['NBRrep','Red','NIR','SWIR2','doy','Redb','NIRb','SWIR2b','doyb']);

var min_NBRdiff = allcoll_diff_select.reduce(ee.Reducer.min(allcoll_diff_select.first().bandNames().size()))
  .rename(['NBRrep','Red','NIR','SWIR2','doy','Redb','NIRb','SWIR2b','doyb']);

print('Image reduced by min NBR and other bands: ', min_NBRdiff);

var min_NBRdiff_pt = min_NBRdiff.clip(geom);

Map.addLayer(min_NBRdiff_pt, {bands: ['SWIR2', 'NIR', 'Red']}, 'min NBR RGB',false);

Map.addLayer(min_NBRdiff_pt, {bands: ['NBRrep']}, 'min NBR',false);



///////////////////////// R-NBR /////////// 
/// Utiliza-se a coleção antes do Min NBR [ltCollection_all] que tem bandas ['NBR','Red','NIR','SWIR2','doy']

// Function to compute NBR change ratio (NBR_t1 - NBR_t2) / (NBR_t2 + 1.001)
function computeNBRChange(previous, current) {
  return previous.subtract(current).divide(current.add(1.001)); // Evita divisão por zero
}


// Function to get annual composite based on max NBR change ratio within the same year, preserving DOY and NBR values
function getAnnualComposite(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = ee.Date.fromYMD(year, 12, 31);
  
  var landsat = ltCollection_all
    .filterDate(start, end)
    .sort('system:time_start');
  
  var list = landsat.toList(landsat.size());
  
  var imagesWithRatio = ee.ImageCollection(ee.List.sequence(1, list.size().subtract(1)).map(function(i) {
    var current = ee.Image(list.get(i));
    var previous = ee.Image(list.get(ee.Number(i).subtract(1))); // Get previous image

    // Ensure pixel-level validity by using masked areas
    var validPrevious = previous.updateMask(previous.select('NBR').mask());
    var validCurrent = current.updateMask(current.select('NBR').mask());
    
    var rNBR = computeNBRChange(validPrevious.select('NBR'), validCurrent.select('NBR')).rename('rNBR');

    // Preserve DOY, NBR previous, and NBR current
    var doyCurrent = current.select('doy').rename('DOY_Current');
    var doyPrevious = previous.select('doy').rename('DOY_Previous');
    var nbrPrevious = validPrevious.select('NBR').rename('NBR_Previous');
    var nbrCurrent = validCurrent.select('NBR').rename('NBR_Current');

    return current.addBands([rNBR, doyCurrent, doyPrevious, nbrPrevious, nbrCurrent]);
  }));
  
  var maxChangeImage = imagesWithRatio.qualityMosaic('rNBR');
  return maxChangeImage;
}

// Generate annual composites for 2019 and 2020
var years = ee.List.sequence(2019, 2020);
var annualComposites = ee.ImageCollection(years.map(getAnnualComposite));

// Get current year (2020) and previous year (2019) composites
var currentYear = 2020;
var previousYear = 2019;

var currentComposite = getAnnualComposite(currentYear);
var previousComposite = getAnnualComposite(previousYear);

print('Annual Composite ' + currentYear ,currentComposite)


// Compute the difference between rNBR of current and previous year
var rNBR_Diff = currentComposite.select('rNBR').subtract(previousComposite.select('rNBR')).rename('rNBR_Diff');


// Apply Haralick homogeneity based on DOY_Current
var doyCurrent = currentComposite.select('DOY_Current').toUint16(); // Convert to integer

var doyCurrentTexture = doyCurrent.glcmTexture({size: 3}).select('DOY_Current_idm'); 
  // das várias bandas criadas a de interesse é : IDM: f5, Inverse Difference Moment; measures the homogeneity

// Add layers to the map

// para ser mais rapida a visualização, criação de .clip(geom)
Map.addLayer(currentComposite.clip(geom), {bands: ['SWIR2', 'NIR', 'Red'], min: 0, max: 0.3}, 'Annual Composite ' + currentYear);
Map.addLayer(previousComposite.clip(geom), {bands: ['SWIR2', 'NIR', 'Red'], min: 0, max: 0.3}, 'Annual Composite ' + previousYear);
Map.addLayer(rNBR_Diff.clip(geom), {min: -1, max: 1, palette: ['blue', 'white', 'red']}, 'rNBR Difference (2020 - 2019)');
  //HOMOGENEIDADE VARIA ENTRE 0 e 1
Map.addLayer(doyCurrentTexture.clip(geom), {min: 0, max: 1, palette: ['black', 'white']}, 'DOY Homogeneity (Haralick)');

