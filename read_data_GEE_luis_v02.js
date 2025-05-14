// -------- FIREATLAS -------- //

// Variaveis a definir: 
// RegionName; currentyear; Satelites (sat1; sat1Name; sat2; sat2Name)

// Outputs:
  // - Tabela com num. de imagens por sensor/satelite
  // - Compósito anual, para metade de PT, com 8 bandas


// ===  REGIÃO a ser processada === 
var regionName = 'Norte';    // definir como 'Norte' ou 'Sul'

// ===  Ano a ser processado === 
var currentyear = 2012;

    // 1984–1993	Landsat 4, Landsat 5
    // 1994–1998	Landsat 5
    // 1999–2002	Landsat 5, Landsat 7
    // 2003–2012	Landsat 5, Landsat 7 (SLC-off desde 2003) *
    // 2013–2021	Landsat 7, Landsat 8
    // 2022–2024	Landsat 8, Landsat 9
    
      // * No caso de 2003–2011, dada a situação do SLC-off, o L7 deve ser 
      // considerado Sat1 e L7 como Sat2, prioritizando imagens do L5
      
      
    // PRE‐DEFINE todas as coleções possiveis:
    var Landsat4  = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2');
    var Landsat5  = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');
    var Landsat7  = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2');
    var Landsat8  = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
    var Landsat9  = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');


// === SATÉLITES a serem usados ===
var sat1 = Landsat8; // ex: Landsat-8 
var sat1Name = 'LANDSAT_8';  // Nome a ser exportado na coleção c/ o nº de imagens

var sat2 = Landsat9;  
var sat2Name = 'LANDSAT_9';




// ===  BANDAS a exportar === 
// Todas as bandas geradas menos Red, NIR, SWIR2 relativas ao pré-fogo
//var bandstoexport = ['Red', 'NIR', 'SWIR2', 'NBRrep', 'NBRrep_b', 'NBRChange', 'doy', 'year', 'doyb', 'yearb', 'doy_homogeneity'];

var bandstoexport = ['Red', 'NIR', 'SWIR2', 'NBRrep', 'NBRChange', 'doy', 'timeDiff', 'doy_homogeneity'];




// ---- Area de interesse ----
// Asset com limites de Portugal continental (com buffer de 1000m), 
//divido em duas partes (Norte e Sul)
var portugal = ee.FeatureCollection('projects/ee-llopes/assets/Portugal_limits_1000_2partes'); 

// Filtrar a coleção para obter apenas a região escolhida
var geom = portugal.filter(ee.Filter.eq('Regiao', regionName));

Map.addLayer(portugal, {color: 'yellow'}, 'ROI');
Map.centerObject(geom, 12);
//Map.addLayer(geometry ,  {color: 'yellow'}, 'POI');




// --------------- PROCESSAMENTO ---------------

// Data de inicio e fim para a recolha de imagens / analise
var start = ee.Date.fromYMD(currentyear, 1, 1).advance(-30, 'day'); // adicionado mes anterior à coleção
var end = ee.Date.fromYMD(currentyear, 12, 31);


// Compilacao de imagens / criacao de colecao para cada landsat

var filteredCollection_sat1 = sat1
  .filter(ee.Filter.bounds(geom))
  .filterMetadata('CLOUD_COVER', 'less_than', 40)  
  .filterDate(start, end).sort('system:time_start');
  
  
var filteredCollection_sat2 = sat2
  .filter(ee.Filter.bounds(geom))
  .filterMetadata('CLOUD_COVER', 'less_than', 40)
  .filterDate(start, end).sort('system:time_start');  
  


// Obter e guardar número de imagens por coleção
print('Número de imagens ' + sat1Name, filteredCollection_sat1.size());
print('Número de imagens ' + sat2Name, filteredCollection_sat2.size());



var Imagensporsatelite = ee.FeatureCollection([
  ee.Feature(null, {
    sensor: sat1Name,
    count: filteredCollection_sat1.size(),
    year: currentyear
  }),
  ee.Feature(null, {
    sensor: sat2Name,
    count: filteredCollection_sat1.size(),
    year: currentyear
  })
]);

Export.table.toDrive({
  collection: Imagensporsatelite,
  folder: 'GEE_EXPORT_ATLAS', 
  description: 'NumImgSatelite_' + currentyear + regionName,
  fileFormat: 'CSV'
});

  

/////////////////////// RENAME Sat1 AND Sat2 BANDS, ADD DOY BAND AND MERGE THEM IN A UNIQUE COLLECTION///////////////////////
  //input: filteredCollection_sat1 AND filteredCollection_sat2
  //output: col_mosaic 
  
  
    // --- RENAME BANDS 
 
    
// Alterar nomes em "Imagensporsatelite" que dá origem à tabela 
// Confirmar se necessárias alteracoes nas coleção, nos nomes/bandas em "Rename bands"
// As duas “combinações” possiveis de bandas:
var combo1 = ['SR_B3','SR_B4','SR_B7','ST_B6','QA_PIXEL','QA_RADSAT'];
//Apenas L8 e L9 têm bandas diferentes para o nosso objetivo   
var combo2 = ['SR_B4','SR_B5','SR_B7','ST_B10','QA_PIXEL','QA_RADSAT'];

// Nomes finais (sempre os mesmos):
var nomesbandas =    ['Red','NIR','SWIR2','Thermal','QA_PIXEL','QA_RADSAT'];    
    
   
   
//Sat1
 // Usa um if/else para apontar bandselect à combo certa
    var bandselect_sat1;
    if (sat1 === Landsat8 ||
        sat1 === Landsat9) {
      bandselect_sat1 = combo2;
    } else {
      bandselect_sat1 = combo1;
    } 
  
  
  // Alteracao dos nomes
  var rename_sat1 = function(img) {
      var select = img.select(bandselect_sat1) 
                      .rename(nomesbandas)
                    .set('system:time_start', img.get('system:time_start'));
      return select;
  };
  
  var Sat1_renamed = filteredCollection_sat1.map(rename_sat1);
  
  //print('Sat1_renamed:', Sat1_renamed);
  //Map.addLayer(Sat1_renamed, {}, 'Sat1_renamed');
  
  
  
  //Sat2
 // Usa um if/else para apontar bandselect à combo certa
    var bandselect_sat2;
    if (sat2 === Landsat8 ||
        sat2 === Landsat9) {
      bandselect_sat2 = combo2;
    } else {
      bandselect_sat2 = combo1;
    } 
  
  
  // Alteracao dos nomes
  var rename_sat2 = function(img) {
      var select = img.select(bandselect_sat2) 
                      .rename(nomesbandas)
                    .set('system:time_start', img.get('system:time_start'));
      return select;
  };
  
  var Sat2_renamed = filteredCollection_sat2.map(rename_sat2);
  
  print('Sat2_renamed:', Sat2_renamed);
  
  Map.addLayer(Sat2_renamed, {}, 'Sat1_renamed');
  
  
 
  // ---  CREATE DATE BAND (DOY) FOR EACH IMAGE IN COLLECTION
  
  var addDate = function(image){
  // Compute Day of Year (DOY)
    var doy = image.date().getRelative('day', 'year'); 
    var doyBand = ee.Image.constant(doy).int16().add(1).rename('doy');
	// Compute Year
    var year = image.date().get('year');
    var yearBand = ee.Image.constant(year).int16().rename('year');
  // Apply based on a band to create a mask 
    var mask = image.select('NIR').mask();// Escolhi 'NIR' como referência para manter consistência
    doyBand = doyBand.updateMask(mask);
    yearBand = yearBand.updateMask(mask);
  // Add both
    return image.addBands([doyBand, yearBand]);
  };
  
  //apply it
  var Sat1_withdate = Sat1_renamed.map(addDate).sort('system:time_start');
  print(sat1Name + ' Collection with doy band', Sat1_withdate);
  
  var Sat2_withdate = Sat2_renamed.map(addDate).sort('system:time_start');
  print(sat2Name + ' Collection with doy band', Sat2_withdate);
 
  Map.addLayer(Sat1_withdate, {}, sat1Name + ' Collection with doy band');
  Map.addLayer(Sat2_withdate, {}, sat2Name + ' Collection with doy band');
 

// ---  MERGE Satelite 1 e Satelite 2
//input: Sat1_withdate e Sat2_withdate
//output: col_mosaic


// Antes do merge, adicionar uma banda de prioridade para, nas situações em que existem imagens para o mesmo dia
// ser dada prioridade as imagens do satelite mais recente
// Anteriormente essa escolha era feita com base na imagem mais recente, independentemente da fonte


// Adicionar a banda priority
// valor 1 para Landsat Sat1 e valor 2 para Landsat Sat

var Sat1_renamed_priority = Sat1_withdate.map(function(img) {
  return img.addBands(ee.Image.constant(1).rename('priority'))
            .set('system:time_start', img.get('system:time_start'));
});

var Sat2_renamed_priority = Sat2_withdate.map(function(img) {
  return img.addBands(ee.Image.constant(2).rename('priority'))
            .set('system:time_start', img.get('system:time_start'));
});


// ---  CREATE MOSAICS WITH IMAGES FROM THE SAME DATE
function mosaicByDate(imcol){
  var imlist = imcol.toList(imcol.size());
  
  var unique_dates = imlist.map(function(im){
    return ee.Image(im).date().format("YYYY-MM-dd");
  }).distinct();
  
  var mosaic_imlist = unique_dates.map(function(d){
    d = ee.Date(d);

    var im = imcol
      .filterDate(d, d.advance(1, "day"))
      .qualityMosaic('priority');  

    return im.set(
      "system:time_start", d.millis(),
      "system:id", d.format("YYYY-MM-dd"));
  });
  
  return ee.ImageCollection(mosaic_imlist);
}


// Uniao das coleções
var filteredCollection_all = Sat2_renamed_priority.merge(Sat1_renamed_priority).sort('system:time_start');

// Cria os mosaicos 
var col_mosaic = mosaicByDate(filteredCollection_all).sort('system:time_start');

// Remover banda 'priority' após o mosaico
var col_mosaic = col_mosaic.map(function(img) {
    return img.select(img.bandNames().remove('priority'));
});
                   
//Map.addLayer(col_mosaic, {}, 'col_mosaic');
print('Mosaics of the same date', col_mosaic);


  
///////////////////////  MASCARA DE QUALIDADE baseada em QA_PIXEL e QA_RADSAT, e thermalBand > 300 ///////////////////////
//input: col_mosaic 
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
  //print('maskedcollection_all',maskedcollection_all);

  
//////////////////////  ADD NBR BAND (NO TEMPORAL DIFFERENCES) ///////////////////////
// input: maskedcollection_all
// output: ltCollection_all_NBR

  // CREATE NBR FUNCTION 
  var nbr = function(img) {
   
      var index = img.normalizedDifference(['NIR', 'SWIR2'])
                    .select([0], ['NBR'])
                    .set('system:time_start', img.get('system:time_start'));
      return img.addBands(index).toFloat() ; 
  };
  
  //Apply it
  var ltCollection_all_NBR =  maskedcollection_all.map(nbr); 
  
//print('Filtered collection with NBR', ltCollection_all_NBR);


/////////////////////// FILTRO DE HAMPEL ///////////////////////

// input: ltCollection_all_NBR (image collection) com banda NBR
// OUTPUT: ltCollection_all_hampel: image collection, bandas ['Red','NIR','SWIR2','doy','NBR','year','Redb','NIRb','SWIR2b','doyb','NBRb','NBRdiff']

	// --- CALCULA NBRmedian  ---------
	// OUTPUT: IC NBRmedian

  //select only NBR bands
// a banda 'NBR' está em  ltCollection_all: porque é necessário calcular a IC NBRdiff
  var nbrbands = ltCollection_all_NBR.select(['NBR']);  // MC: não a dif de NBR , pois é a banda 'NBR' de NBRdiff, não e banda 'NBRdiff' // LL substitui NBRdiff por ltCollection_all

// MC: calcula a mediana de 3 sucessivas imagens? Sim
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

	// ---  substitui outliers em nbrbands por NBRmedian -----
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
    //print('replaced',replaced);
  //Map.addLayer(replaced, {}, 'replaced');

  // Combina as bandas processadas com a coleção original
// Nesta abordagem, para cada imagem original é adicionado o NBRmedian, NBRmedian_diff e NBRrep
var ltCollection_all_hampel = ltCollection_all_NBR.map(function(img) {
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
  Map.addLayer(ltCollection_all_hampel, {}, 'ltCollection_all_hampel');



	///////////  R_NBR ////////////
// input: ltCollection_all_hampel
// output: ltCollection_upd_RNBR

//CALCULAR o ÍNDICE R_NBR que depois é usado para o REDUCER  
//Utiliza a banda NBRrep que é a banda NBR atualizada apos o filtro de hampel com 

// Bandas "atuais"
var bands = ['Red','NIR','SWIR2','doy','year','NBRrep'];
// Bandas _b - before -  que irão conter o "carry-forward" (valores anteriores)
var bands_b = ['Redb','NIRb','SWIR2b','doyb','yearb','NBRrep_b'];

// Ordena sua coleção final (pós-Hampel) por data e converte para lista
var sorted = ltCollection_all_hampel.select(bands).sort('system:time_start');
var list   = sorted.toList(sorted.size());

// Imagem dummy inicial: sem valores válidos (0), para a primeira iteração
var dummyImage = ee.Image.constant([0,0,0,0,0,0])
  .rename(bands)
  .toFloat();


//////////  FUNÇÃO DE ITERAÇÃO

// Nesta função, `accumImg` é a imagem que contém o "carry-forward"
// de todas as iterações passadas. Em cada iteração:
//   - Gera-se 'combined' = (imagem atual + bands_b do accumImg).
//   - Atualiza-se 'accumImg' com os valores válidos da imagem atual
//     (assim, se a próxima imagem estiver mascarada em certo pixel,
//     continua-se usando o valor antigo).


var iterateFunction = function(i, state) {
  i     = ee.Number(i);
  state = ee.Dictionary(state);
  
  // "accumImg": última versão do carry-forward
  var accumImg = ee.Image(state.get('accumImg'));
  // "results": lista de imagens combinadas já geradas
  var results  = ee.List(state.get('results'));
  
  // Pega a imagem atual (sem bandas_b), forçando float
  var current = ee.Image(list.get(i)).toFloat();
  
  // imagem combinada: bandas atuais = current / bandas_b = accumImg (renomeado)
  var combined = current.addBands(accumImg.rename(bands_b))
    .copyProperties(current, current.propertyNames());
  
  // carry-forward para cada pixel:
  //       - Se current estiver válido, sobrescreve accumImg
  //       - Se current estiver mascarado, mantém o valor anterior
  var updatedAccum = current.unmask(accumImg);
  
  // Atualiza a lista de resultados e retorna o novo estado
  var newResults = results.add(combined);
  
  return ee.Dictionary({
    'accumImg': updatedAccum, 
    'results': newResults
  });
};


//////////  EXECUTA A ITERAÇÃO 

// Cria uma lista de índices [0..n-1]
var indices = ee.List.sequence(0, list.size().subtract(1));

// Define o estado inicial: carry-forward = dummy, resultados = []
var initialState = ee.Dictionary({
  'accumImg': dummyImage,
  'results': ee.List([])
});

// Executa a iteração 
var finalState = indices.iterate(iterateFunction, initialState);
var finalDict  = ee.Dictionary(finalState);

// Extrai a lista de imagens combinadas
var resultsList = ee.List(finalDict.get('results'));

// Converte para ImageCollection 
var ltCollection_upd_RNBR = ee.ImageCollection(resultsList);


//Map.addLayer(ltCollection_upd_RNBR, {bands: ['SWIR2b','NIRb','Redb']}, 'Valores_b (carry-forward)');
//print('withBeforeColl', ltCollection_upd_RNBR);




// ---  Function to compute NBR change ratio  ------
	
	// Formula/Ratio = (dNBR)/(NBR + 1.001)
    // (NBR_Prefire - NBR_Postfire) / (NBR_Postfire + 1.001)
	
// input: ltCollection_upd_RNBR
// output: ltCollection_upd_NBRchange


// Manually define leap years from 1984 to 2024
var leapYears = ee.List([
  1988, 1992, 1996, 2000, 2004, 2008,
  2012, 2016, 2020, 2024
]);
//print('Leap years:', leapYears);


// Determine if current year is a leap year
var isLeap = leapYears.contains(currentyear);
var daysInYear = ee.Number(ee.Algorithms.If(isLeap, 366, 365));
print(daysInYear)


var addNBRChange = function(image) {
  // Aplica a máscara apenas à imagem atual
  var yearMask = image.select('year').eq(currentyear);

  // dNBR = (NBR_Prefire - NBR_Postfire)
  var dNBR = image.select('NBRrep_b').subtract(image.select('NBRrep'));

 // Diferença base (caso o ano seja o mesmo)
    var timeDiff = image.select('doy').subtract(image.select('doyb'));
    
    // Corrige para casos em que os anos são diferentes
    timeDiff = timeDiff.where(
      image.select('year').neq(image.select('yearb')),
      image.select('doy').subtract(image.select('doyb').subtract(daysInYear))
    ).rename('timeDiff');
    
  // Calcula NBRChange apenas onde year == currentyear
  var NBRChange = dNBR
    .divide(image.select('NBRrep').add(1.001))
    .updateMask(yearMask)
    .rename('NBRChange');

  return image.addBands([timeDiff, NBRChange]);
};

// Map the function over the collection to get the final result
var ltCollection_upd_NBRchange = ltCollection_upd_RNBR.map(addNBRChange);

Map.addLayer(ltCollection_upd_NBRchange, {bands: ['SWIR2','NIR','Red']}, 'ltCollection_upd_NBRchange');



////////// Multiply All Pixels by 10,000 and Convert to Integer (temporal bands not multiplied)

// List of bands to exclude from scaling
var bandsToKeep = ['doy', 'year','doyb', 'yearb', 'timeDiff'];

// Function to scale only selected bands and keep others as integer
var scaleAndConvert = function(img) {
    var bandsToScale = img.bandNames().removeAll(bandsToKeep); // Get all bands except those in 'bandsToKeep'
    
    var scaledBands = img.select(bandsToScale)  // Select bands to scale
                         .multiply(10000)       // Scale values
                         .toInt16();            // Convert to int16

    var unscaledBands = img.select(bandsToKeep)  // Select bands to keep unchanged
                            .toInt16();          // Convert to integer format

    return scaledBands.addBands(unscaledBands)   // Merge both sets of bands
                      .copyProperties(img, img.propertyNames());  // Keep metadata
};

// Apply the function to all images in the collection
var NBRChange_scaled = ltCollection_upd_NBRchange.map(scaleAndConvert);
Map.addLayer(NBRChange_scaled, {bands: ['SWIR2','NIR','Red']}, 'NBRChange_scaled');


// Print the scaled collection to check
print('NBRChange_scaled', NBRChange_scaled);


///////////////////// REDUCER PARA OBTER COMPOSITO ANUAL //////////////////////////

// REDUCE BY MAX_NBRCHANGE 
// "The qualityMosaic() function selects the image (per-pixel) with the HIGHEST quality-band-score to contribute to the resulting mosaic". That is, the values are taken from the image where NBRChange is maximum.

var MAX_NBRCHANGE = NBRChange_scaled.qualityMosaic('NBRChange');
Map.addLayer(MAX_NBRCHANGE, {bands: ['SWIR2','NIR','Red'], min: -10000, max: 10000}, 'MAX_NBRCHANGE ');


// Print and check the result
//print('MAX_NBRCHANGE composite:', MAX_NBRCHANGE);
//Map.addLayer(MAX_NBRCHANGE, {bands: ['SWIR2','NIR','Red'], min: -10000, max: 10000}, 'MAX_NBRCHANGE ');


///////////////////// HOMOGENEITY //////////////////////////
// input: MAX_NBRCHANGE
// output: AnualComposite

// Apply Haralick homogeneity based on DOY_Current
var doyCurrent = MAX_NBRCHANGE.select('doy').toInt16(); // Convert to integer as Image.glcmTexture: Only 32-bit or smaller integer types are currently supported.

//print('Annual Composite doyCurrent ' + currentyear ,doyCurrent);

var doyCurrentTexture = doyCurrent.glcmTexture({size: 3}).select('doy_idm')
              .rename('doy_homogeneity')
              .multiply(10000)
              .toInt16(); 
  // das várias bandas criadas a de interesse é : IDM: Inverse Difference Moment; measures the homogeneity

//Map.addLayer(doyCurrentTexture.select('doy_homogeneity'), {}, 'doy_homogeneityteste ' + currentyear + 'NBRChange');


// Add the selected band to MAX_NBRCHANGE (currentComposite)
var AnualComposite = MAX_NBRCHANGE.addBands(doyCurrentTexture);

// Seleciona apenas as bandas da imagem final a exportar
var AnualComposite = AnualComposite.select(bandstoexport);


print('Annual Composite ' + currentyear ,AnualComposite);
Map.addLayer(AnualComposite, {bands: ['SWIR2','NIR','Red'], min: 0, max: 10000}, 'AnualComposite ' + currentyear );



//////      EXPORT      //////

//Export MULTIBAND RASTER
 Export.image.toDrive({
   image:AnualComposite,
   description: 'CompositoAnual_'+currentyear+regionName,
   folder: 'GEE_EXPORT_ATLAS',
   scale: 30,
   region: geom,
   maxPixels: 1e13,
   fileFormat: 'GeoTIFF',
   crs: 'EPSG:4326' // -> crs WGS84 
 });
 
 
