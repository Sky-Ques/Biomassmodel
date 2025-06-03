
// ─── (A) Denmark boundary ───
var denmark = ee.FeatureCollection("FAO/GAUL/2015/level0")
              .filter(ee.Filter.eq('ADM0_NAME', 'Denmark'));
Map.addLayer(denmark, {color:'blue', fillColor:'00000000'}, 'Denmark (GADM)');

// ─── (B) Load & Composite Landsat 8 SR for 2014 ───
// Note: this comes right after  Denmark outline layer.

var landsatSR = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterBounds(denmark)                  // Clip to Denmark
  .filterDate("2014-01-01", "2014-12-31") // All of 2014 (or adjust to leaf-on)
  .select(["SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","QA_PIXEL"]);

print("Landsat 8 SR Collection:", landsatSR);

// Cloud-masking function using QA_PIXEL bit 3
function maskClouds(image) {
  var qa = image.select("QA_PIXEL");
  var mask = qa.bitwiseAnd(1 << 3).eq(0);
  return image.updateMask(mask);
}

// Apply mask & build median composite
var landsatMasked = landsatSR.map(maskClouds);
print("First Landsat Image (masked):", landsatMasked.first());

var landsatImage = landsatMasked
  .median()
  .select(["SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7"])
  .clip(denmark);

print("Final Landsat Composite (2014):", landsatImage);

// Visualization (true-color)
var visParams = {
  bands: ["SR_B4","SR_B3","SR_B2"],
  min: 7000,
  max: 20000,
  gamma: 1.4
};

Map.addLayer(landsatImage, visParams, "L8 2014 Cloud-Free");

// ─── C) 2014 Summer NDVI ───

// 1) Build an NDVI collection for June–Aug 2014
var ndvi2014 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(denmark)
  .filterDate('2014-06-01','2014-08-31')
  .map(maskClouds)
  .map(function(img) {
    return img
      .normalizedDifference(['SR_B5','SR_B4'])
      .rename('NDVI');
  })
  .median()
  .clip(denmark);

// 2) Visualize it
Map.addLayer(
  ndvi2014,
  {min: 0.0, max: 1.0, palette: [
    'FFFFFF', // non‐vegetation
    'CE7E45',
    'DF923D',
    'F1B555',
    'FCD163',
    '99B718',
    '74A901',
    '66A000',
    '529400',
    '3E8601',
    '207401',
    '056201',
    '004C00'  // dense vegetation
  ]},
  'C) NDVI Jun–Aug 2014'
);



// ============================
// Step 1: Load and Process Data
// ============================

// Step 1.1: Load NFI dataset
var nfiData = ee.FeatureCollection("projects/skyques/assets/NFI_2014_LatLon");

// Fix missing geometries by defining points using longitude and latitude columns
var nfiDataFixed = nfiData.map(function(feature) {
  var lon  = ee.Number(feature.get("longitude"));
  var lat  = ee.Number(feature.get("latitude"));
  var geom = ee.Geometry.Point([lon, lat]); 
  return feature.setGeometry(geom);
});

// Print to check if it now has valid geometries
print("Fixed NFI Data:", nfiDataFixed);

// Step 1.2: Load Landsat 8 Surface Reflectance Data (2014) with Cloud Masking
var landsatSR = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterBounds(nfiDataFixed)
  
  .filterDate("2014-01-01", "2014-12-31")
  .select(["SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","QA_PIXEL"]);

// Function to apply cloud masking
function maskClouds(image) {
  var qa    = image.select("QA_PIXEL");
  var mask  = qa.bitwiseAnd(1 << 3).eq(0); // clear‐sky bit
  return image.updateMask(mask);
}




// Apply cloud masking and create a median composite image
var composite = landsatSR
  .map(maskClouds)
  .median()
  .select(["SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7"]);




// ——————————————
// Resampling option
// ——————————————

// Option A: “Nearest‑neighbor” on the native 30 m grid
var landsatNN = composite;  
// (or, if want to be explicit:)
 // var landsatNN = composite
 //   .reproject({
 //     crs: composite.projection(),
 //     scale: 30
 //   });

// Choose nearest‐neighbor for the rest of the pipeline:
var landsatImage = landsatNN;

// ——————————————
// Optional 30 m focal‑mean smoothing
// ——————————————
var originalBands = landsatImage.bandNames();
landsatImage = landsatImage
  .reduceNeighborhood({
    reducer: ee.Reducer.mean(),
    kernel: ee.Kernel.square(1, 'pixels')
  })
  .rename(originalBands);

print("Final Landsat Composite (NN + smoothed):", landsatImage);


// ============================ Step 1.3: Compute Vegetation Indices ============================ //

var landsatWithIndices = landsatImage
  .addBands(landsatImage.normalizedDifference(["SR_B5", "SR_B4"]).rename("NDVI")) // NDVI
  .addBands(landsatImage.expression( // EVI Calculation
    "2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))", {
      "NIR": landsatImage.select("SR_B5"),
      "RED": landsatImage.select("SR_B4"),
      "BLUE": landsatImage.select("SR_B2")
    }).rename("EVI"))
  .addBands(landsatImage.normalizedDifference(["SR_B5", "SR_B7"]).rename("NBR")) // NBR
  .addBands(landsatImage.expression( // SAVI Calculation
    "((NIR - RED) / (NIR + RED + L)) * (1 + L)", {
      "NIR": landsatImage.select("SR_B5"),
      "RED": landsatImage.select("SR_B4"),
      "L": 0.5
    }).rename("SAVI"))
  .addBands(landsatImage.normalizedDifference(["SR_B3", "SR_B5"]).rename("NDWI")) // NDWI
  .addBands(landsatImage.normalizedDifference(["SR_B5", "SR_B6"]).rename("NDMI")) // NDMI
  .addBands(landsatImage.select("SR_B5").divide(landsatImage.select("SR_B3")).rename("GCI")) // GCI

// Update reference to use new Landsat image with indices
// ——— Step 1.3.1: Mask out everything with NDVI ≤ 0 to keep only vegetated pixels ———---------------

var landsatProcessed = landsatWithIndices.updateMask(
  landsatWithIndices.select("NDVI").gt(0)
);
print("→ landsatProcessed (NDVI > 0):", landsatProcessed);


// ————————————————————————————————
// Step 1.3.2: Ingest DEM & derive terrain metrics
// ————————————————————————————————

// 1. Clip & load the 30 m SRTM DEM over study area
var dem = ee.Image("USGS/SRTMGL1_003")
  .clip(nfiDataFixed.geometry().bounds());

// 2. Compute elevation, slope and aspect
var terrain   = ee.Terrain.products(dem);
var elevation = terrain.select("elevation").rename("elev");
var slope     = terrain.select("slope")    .rename("slope");
var aspect    = terrain.select("aspect")   .rename("aspect");

// 3. Convert aspect → northness & eastness
var rad       = aspect.multiply(Math.PI/180);
var northness = rad.cos().rename("northness");
var eastness  = rad.sin().rename("eastness");

// 4. Reproject everything onto  Landsat composite’s 30 m grid
var P = landsatImage.projection();
elevation = elevation.reproject(P);
slope     = slope    .reproject(P);
northness = northness.reproject(P);
eastness  = eastness .reproject(P);

// 5. Add all terrain bands into  processed stack
landsatProcessed = landsatProcessed
  .addBands([elevation, slope, northness, eastness]);

print("→ landsatProcessed (+ DEM & terrain metrics):", landsatProcessed);

//  CHECK: make sure no NFI points dropped out
var nPts = landsatProcessed.sampleRegions({
    collection: nfiDataFixed,
    properties: ["BMag_ha"],
    scale: 30,
    geometries: false
  })
  .filter(ee.Filter.gt("BMag_ha", 0))
  .size();
print("Points surviving with DEM bands:", nPts);


// ——————————————
// Step 2, Vegetation Indices 
// ——————————————

// Step 2.1: Compute Vegetation Indices (NDVI, EVI, NBR, SAVI)
var landsatWithNDVI = landsatImage
  .addBands(landsatImage.normalizedDifference(["SR_B5", "SR_B4"]).rename("NDVI")) // NDVI
  .addBands(landsatImage.expression( // EVI Calculation
    "2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))", {
      "NIR": landsatImage.select("SR_B5"),
      "RED": landsatImage.select("SR_B4"),
      "BLUE": landsatImage.select("SR_B2")
    }).rename("EVI"))
  .addBands(landsatImage.normalizedDifference(["SR_B5", "SR_B7"]).rename("NBR")) // NBR
  .addBands(landsatImage.expression( // SAVI Calculation
    "((NIR - RED) / (NIR + RED + L)) * (1 + L)", {
      "NIR": landsatImage.select("SR_B5"),
      "RED": landsatImage.select("SR_B4"),
      "L": 0.5
    }).rename("SAVI"));

// ============================ Step 3: Sample Training Data (Log-Transformed) ============================ //

// Step 3.1: Sample Landsat at NFI Points and Apply Log-Transformation
var trainingData = landsatProcessed.sampleRegions({
  collection: nfiDataFixed,
  properties: ["BMag_ha"],
  scale: 30,
  geometries: true
}).filter(ee.Filter.gt("BMag_ha", 0)) // Remove missing and zero values
  .map(function(feature) {
    var biomass = ee.Number(feature.get("BMag_ha"));
    var logBiomass = biomass.log(); // Log-transform
    return feature.set("log_BMag_ha", logBiomass);
  });

// Debug: Check first sample after transformation
print("First Training Sample with Log Biomass:", trainingData.first());

// ————————————————
//  CHECK: Inspect training set
// ————————————————
print("TrainingData count:", trainingData.size());
print("First few training samples:", trainingData.limit(5));

// ============================ Step 4: Train Log-Transformed Model ============================ //

// Step 4.2: Apply Grid-Based Stratified Sampling
var gridSize = 10000;
var grid = ee.Image.random().multiply(1000000).toInt();

var stratifiedNFI = trainingData.map(function(feature) {
  var gridID = grid.sample(feature.geometry(), 30).first().get("constant");
  return feature.set("gridID", gridID);
});


// Step 4.3: Split 80% training, 20% validation
var trainingNFI = stratifiedNFI.randomColumn("rand").filter(ee.Filter.lt("rand", 0.8));
var validationNFI = stratifiedNFI.randomColumn("rand").filter(ee.Filter.gte("rand", 0.8));


// (after stratified split code)
print('Number of training points (80 %):', trainingNFI.size());
print('Number of validation points (20 %):', validationNFI.size());



// Step 4.4: Train the Model Using Log-Transformed Biomass
var predictorBands = ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "NDVI", "EVI", "NBR", "elev", "northness","eastness"];

// Step 4.5: Tune Random Forest hyperparameters
var biomassModel = ee.Classifier.smileRandomForest(
    /* numberOfTrees:       */ 800        //800 best,de andre giver worse
    /* variablesPerSplit:   */ //4,         //4 best
    /* minLeafPopulation:   */ //2,          // 2 giver bedste    har prøvet: 5
    /* bagFraction:         */ //0.8,        //0.8
    /* maxNodes:            */ //200      // ik null men200, bedre ved at skrue ned.
)
.setOutputMode("REGRESSION")
.train({
  features: trainingNFI,
  classProperty: "log_BMag_ha",
  inputProperties: predictorBands
});

// Optional: see which predictors mattered most
print("RF explain:", biomassModel.explain());

var explain = biomassModel.explain();
print("Variable importances:", explain.get("importance"));

// ============================ Step 4.6: Predict Biomass ============================ //


// 1) Klassificér og få log-skalaen (ingen bounds-klip her)
var logBiomassPrediction = landsatProcessed
  .classify(biomassModel)
  .rename("log_predicted_Biomass");

// 2) Konverter tilbage til ton/ha
var biomassPrediction = logBiomassPrediction
  .exp()
  .rename("Predicted_Biomass");

// 3) Clip til Danmark-polygonen
var biomassDenmark = biomassPrediction.clip(denmark);


// ============================
// D) Full-Denmark Predicted Biomass (land only via NDWI + Viridis palette)
// ============================ //

// 1) Predict & back-transform
var logBiomass = landsatProcessed
  .classify(biomassModel)
  .rename('log_predicted_Biomass');

var biomass = logBiomass
  .exp()
  .rename('Predicted_Biomass');

// 2) Clip til Danmark
var biomassDK = biomass.clip(denmark);

// 3) Beregn NDWI og maskér vand
var ndwi = landsatImage
  .normalizedDifference(['SR_B3','SR_B5'])
  .rename('NDWI');
var landMask = ndwi.lt(0);
var biomassLand = biomassDK.updateMask(landMask);

// 4) (re)definer visualization-parametre (Viridis uden gamma)
var biomassVis = {
  min: 0,
  max: 350,
  palette: [
    '440154','3b528b','21918c','5cc863','fde725'
  ]
};

// 5) Tegn laget og center
Map.addLayer(biomassLand, biomassVis, 'D) Biomasse (0–350 t/ha, Viridis)');
Map.centerObject(denmark, 6);




// ============================ Step 5: Compute R² and RMSE ============================ //

// Step 5.1: Extract actual and predicted values
var validationPredictions = logBiomassPrediction.sampleRegions({
  collection: validationNFI,
  properties: ["BMag_ha"],
  scale: 30,
  geometries: true
}).map(function(feature) {
    var logPredicted = ee.Number(feature.get("log_predicted_Biomass"));
    var predictedBiomass = logPredicted.exp();
    return feature.set("predicted_BMag_ha", predictedBiomass);
  });

// Compute R²
var regressionStats = validationPredictions
    .reduceColumns(ee.Reducer.pearsonsCorrelation(), ["BMag_ha", "predicted_BMag_ha"]);

var rSquared = ee.Number(regressionStats.get("correlation")).pow(2);
print("R² (Coefficient of Determination):", rSquared);

// Compute squared errors correctly
var squaredErrors = validationPredictions.map(function(feature) {
    var diff = ee.Number(feature.get("BMag_ha")).subtract(feature.get("predicted_BMag_ha"));
    return feature.set("squaredError", diff.pow(2));
});

// Reduce to mean squared error
var mse = squaredErrors.aggregate_mean("squaredError");
var rmse = mse.sqrt();
print("RMSE (Root Mean Square Error):", rmse);


// ============================ Step 6: Visualization ============================ //   har // på visualisering af residualer

// Add Predicted Biomass Layer
var validationVis = {min: 0, max: 350, palette: ["black", "brown", "darkgreen", "lightgreen", "yellow"]};
Map.addLayer(biomassPrediction, validationVis, "Predicted Biomass");

// Compute residuals (prediction errors)
var residuals = validationPredictions.map(function(feature) {
  var actual = ee.Number(feature.get("BMag_ha"));
  var predicted = ee.Number(feature.get("predicted_BMag_ha"));
  return feature.set("residual", actual.subtract(predicted).abs());
});

// Find locations with highest errors
//var worstErrors = residuals.sort("residual", false).limit(10);
//Map.addLayer(residuals, {color: "red"}, "Residual Errors");
//Map.addLayer(worstErrors, {color: "black"}, "Worst Errors");

// Step 6.1: Generate Scatterplot
var chart = ui.Chart.feature.byFeature(validationPredictions, "BMag_ha", "predicted_BMag_ha")
  .setChartType("ScatterChart")
  .setOptions({
    title: "Actual vs Predicted Biomass",
    hAxis: {title: "Actual Biomass (BMag_ha)"},
    vAxis: {title: "Predicted Biomass"},
    pointSize: 3,
    trendlines: { 0: {color: "red"} }
  });

print(chart);

//-------------------------------------------- Analyze Model Errors  Inspect residuals (errors)


// Step 1.1: Compute residual stats
var residualStats = residuals.aggregate_stats("residual");
print("Residual Statistics:", residualStats);

// Visualize residuals with a color ramp
//var residualVis = {
//  min: 0,
//  max: 150, // adjust as needed based on stats
//  palette: ["green", "yellow", "red"]
//};
//Map.addLayer(residuals, residualVis, "Residuals (Color by Error)");

// Step 1.2: Overlay worst points
//Map.addLayer(worstErrors.style({color: "black", pointSize: 5}), {}, "Worst Residuals");

// Optional: Zoom to these points manually or use high-res background imagery




// ============================ Step 7: Single Pixel Biomass Tracking (Full Time Series - Handling Landsat 5 Bands) ============================ //

// 7.1: Define Pixel of Interest and Time Range
var targetPixel = ee.Geometry.Point([9.226201, 56.822027]);     //    Trend skov: Pixel 1: 9.226201, 56.822027   Pixel 2: 9.238934, 56.812888    
var startYear = 1985;                                                           
var endYear   = 2024; // Set the full time range  

Map.centerObject(targetPixel, 15);
Map.addLayer(targetPixel, {color: 'red'}, 'Target Pixel');

print('Target Pixel:', targetPixel);
print('Start Year:', startYear);
print('End Year:', endYear);



// ============================
// Step 7.2: Function to process Landsat imagery for a single year
// ============================


// new: add a 4th flag `returnImage`
var processLandsatYear = function(landsatSR_year_collection, currentYear, targetPixel, returnImage) {


  // ——————————————————————————————
  // 1) Build two medians (5/7 vs 8) then pick server-side
  var cloudFree = ee.ImageCollection(landsatSR_year_collection).map(maskClouds);
  var base      = ee.Image(cloudFree.median());

  // Landsat 5/7 median (bands SR_B1–B5 + SR_B7)
  var median57 = base
    .select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'])
    .reduceNeighborhood({
      reducer: ee.Reducer.mean(),
      kernel: ee.Kernel.square(1)
    })
    .rename(['BLUE','GREEN','RED','NIR','SWIR1','SWIR2']);

  // Landsat 8 median (bands SR_B2–SR_B7)
  var median8 = base
    .select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'])
    .reduceNeighborhood({
      reducer: ee.Reducer.mean(),
      kernel: ee.Kernel.square(1)
    })
    .rename(['BLUE','GREEN','RED','NIR','SWIR1','SWIR2']);

  // Choose the correct one on the server
  var medianImage = ee.Image(ee.Algorithms.If(
    currentYear.lt(2013),
    median57,
    median8
  ));
  
  
  // ——————————————————————————————
  // 2) Compute indices on that median
  var withIdx = medianImage
    .addBands(medianImage.normalizedDifference(['NIR','RED'])
                        .rename('NDVI'))
    .addBands(medianImage.expression(
      '2.5*((NIR-RED)/(NIR+6*RED-7.5*BLUE+1))', {
        NIR: medianImage.select('NIR'),
        RED: medianImage.select('RED'),
        BLUE: medianImage.select('BLUE')
      }).rename('EVI'))
    .addBands(medianImage.normalizedDifference(['NIR','SWIR2'])
                        .rename('NBR'))
    .addBands(medianImage.expression(
      '((NIR-RED)/(NIR+RED+0.5))*(1+0.5)', {
        NIR: medianImage.select('NIR'),
        RED: medianImage.select('RED')
      }).rename('SAVI'))
    .addBands(medianImage.normalizedDifference(['GREEN','NIR'])
                        .rename('NDWI'))
    .addBands(medianImage.normalizedDifference(['NIR','SWIR1'])
                        .rename('NDMI'))
    .addBands(medianImage.select('NIR')
                        .divide(medianImage.select('GREEN'))
                        .rename('GCI'));
  

  // 3) DEM & terrain metrics (no tiny‐buffer clip)
  var dem = ee.Image('USGS/SRTMGL1_003');
  var terrain = ee.Terrain.products(dem);
  var elev = terrain.select('elevation').rename('elev')
               .reproject(withIdx.projection());
  var aspectRad = terrain.select('aspect').rename('aspect')
                   .reproject(withIdx.projection())
                   .multiply(Math.PI/180);
  var northness = aspectRad.cos().rename('northness');
  var eastness  = aspectRad.sin().rename('eastness');
  
  var fullStack = withIdx.addBands([elev, northness, eastness]);
  
  // 4) Remap generic bands → the RF’s SR_B* input names
  var generic = [
    'BLUE','GREEN','RED','NIR','SWIR1','SWIR2',
    'NDVI','EVI','NBR','elev','northness','eastness'
  ];
  var rfBands = [
    'SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7',
    'NDVI','EVI','NBR','elev','northness','eastness'
  ];
  var ready = fullStack.select(generic, rfBands);
  
  // 5) Classify (log‐biomass) → exp → pull pixel
  var logPred = ready.classify(biomassModel)
                     .rename('log_predicted_Biomass');
  var predImg = logPred.exp().rename('Predicted_Biomass');

  // — if caller wants the image, return it now — 
  if (returnImage === true) {
    return predImg;
  }


//Sample the pixel and guard against empty result instead of the section commented out
  
  //return predImg.reduceRegion({
  //reducer: ee.Reducer.first(),
  //geometry: targetPixel,
  //scale: 30
//}).get('Predicted_Biomass');


//-----------
  
  //Old pixel thing
  //var samples = predImg.sample({
    //region: targetPixel.buffer(90),
    //scale: 30,
    //numPixels: 1
  //});

  // If there is at least one feature, take its 'Predicted_Biomass';
  // otherwise return null.
  //var pb = ee.Algorithms.If(
    //samples.size().gt(0),
    //samples.first().get('Predicted_Biomass'),
    //null
  //);


    // new pixel thing Pull the pixel directly by reduceRegion (no chance of random sampling miss).
  //var dict = predImg.reduceRegion({
    //reducer: ee.Reducer.first(),
    //geometry: targetPixel,
    //scale: 30
  //});
  //return ee.Algorithms.If(
    //dict.contains('Predicted_Biomass'),
    //dict.get('Predicted_Biomass'),
    //null
  //);

    // NEW: 3×3 window mean
  var samples = predImg.sample({
    region: targetPixel.buffer(45),  // 45 m radius ≈ 3×3 pixels   3×3 window (nine 30 m pixels) around  target point
    scale: 30,
    numPixels: 9
  });
  var pb = ee.Number(samples.aggregate_mean('Predicted_Biomass'));

  return pb;
//------------------
};


// ——— MINI TEST FOR 2012–2015 ———
var testYears = [2012, 2013, 2014, 2015, 2023, 2024];
ee.List(testYears).evaluate(function(years) {
  years.forEach(function(y) {
    var yr = ee.Number(y);


    // Pick the same collections used later:
    var coll = ee.ImageCollection(
      ee.Algorithms.If(
        yr.gte(2013), col13_on,
      ee.Algorithms.If(
        yr.gte(1999), col99_12,
                       col85_98
      ))
    )
    .filterDate(ee.Date.fromYMD(yr,5,15), ee.Date.fromYMD(yr,9,15))
    .filterBounds(targetPixel);

    // Count how many make it through mask:
    var sceneCount = coll.map(maskClouds).size();

    // Run per-year processor:
    var biomass = processLandsatYear(coll, yr, targetPixel);
    

    
        // inside  ee.List(testYears).evaluate(...) loop
    var sceneCount = coll
      .map(maskClouds)
      .size();
    print('Year', y, '→ scenes:', sceneCount, '→ biomass:', biomass);

    
  });
});



// ============================
// Step 7.3: Process all years & extract the pixel values
// ============================

// Pre-define the three server-side ImageCollection objects
var col85_98 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');
var col99_12 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2');
var col13_on = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');



var biomassTimeSeriesPixel = ee.FeatureCollection(
  ee.List.sequence(startYear, endYear).map(function(y) {
    y = ee.Number(y);
    var yearInt    = y.toInt();
    var startDate = ee.Date.fromYMD(yearInt, 5, 15);
    var endDate    = ee.Date.fromYMD(yearInt, 9, 15);

    // Pick the right collection object:
    var collYear = ee.ImageCollection(
      ee.Algorithms.If(
        y.gte(2013), col13_on,
      ee.Algorithms.If(
        y.gte(1999), col99_12,
                      col85_98
      )))
      .filterDate(startDate, endDate)
      .filterBounds(targetPixel);

    print('Year:', yearInt, 'Image Collection Size:', collYear.size());

    // — build a Feature carrying both mean & stdDev — 
    return ee.Feature(
      ee.Algorithms.If(
        collYear.size().gt(0),
        (function() {
          // 1) get the full image
          var img = processLandsatYear(collYear, y, targetPixel, true);

          // 2) reduceRegion for mean + stdDev over 3×3 window
          var stats = img.reduceRegion({
            reducer: ee.Reducer
                      .mean()
                      .combine({
                        reducer2: ee.Reducer.stdDev(),
                        sharedInputs: true
                      }),
            geometry: targetPixel.buffer(45),
            scale: 30,
            bestEffort: true
          });

          // 3) return feature with both values
          return ee.Feature(null, {
            year:            y,
            biomass_mean:    stats.get('Predicted_Biomass_mean'),
            biomass_stdDev:  stats.get('Predicted_Biomass_stdDev')
          });
        })(),
        // no scenes - nulls
        ee.Feature(null, {
          year:           y,
          biomass_mean:   null,
          biomass_stdDev: null
        })
      )
    );

  })
);

print('biomassTimeSeriesPixel:', biomassTimeSeriesPixel);

var missing = biomassTimeSeriesPixel.filter(ee.Filter.notNull(['predicted_biomass']));
var hasNulls = biomassTimeSeriesPixel.filter(ee.Filter.eq('predicted_biomass', null));
print('Years with value:', missing);
print('Years without value:', hasNulls);

// Export to Drive
Export.table.toDrive({
  collection: biomassTimeSeriesPixel,
  description: 'single_pixel_biomass_1985_2024',
  folder: 'EarthEngineExports',
  fileNamePrefix: 'single_pixel',
  fileFormat: 'CSV'
});



print('Exporting single pixel biomass time series table to Google Drive');

//Should be containing the year and the predicted biomass for the target pixel for each year
//        (with potential null values for years with no suitable imagery).





