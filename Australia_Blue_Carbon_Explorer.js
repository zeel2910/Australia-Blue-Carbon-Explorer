
 * Australia Blue-Carbon Explorer

/*  CONFIG  */
var ASSETS = {
  mangroves:   'users/zeelchavda911/Mangroves',
  seagrasses:  'users/zeelchavda911/Seagrasses',
  saltmarshes: 'users/zeelchavda911/Saltmarshes'
};

var YEARS = (function(){ var arr=[]; for (var y=2001; y<=2024; y+=5) arr.push(String(y)); if (arr[arr.length-1] !== '2024') arr.push('2024'); return arr; })();

var PALETTES = {
  ndvi:  ['#d7191c','#fdae61','#ffffbf','#abdda4','#2b83ba'],
  dndvi: ['#d73027','#fefefe','#1a9850']
};

var THEME = {
  cardBg:  'rgba(255,255,255,0.97)',
  cardHdr: '#f5f7fa',
  text:    '#1f2937',
  subText: '#6b7280',
  border:  '#e5e7eb',
  accent:  '#334155'
};

var STYLES = {
  mangrove:  {color:'1b9e77', fillColor:'1b9e7740', width:1},
  seagrass:  {color:'377eb8', fillColor:'377eb840', width:1},
  saltmarsh: {color:'e6ab02', fillColor:'e6ab0240', width:1}
};

var STATE_CANDIDATES = [
  'ste_name21','STE_NAME21','STATE','STATE_NAME','State_Name','state','STATE_NAME21'
];

var GAUL_NAME_MAP = { // tolerant mapping
  'NSW':'New South Wales','New South Wales':'New South Wales',
  'QLD':'Queensland','Queensland':'Queensland',
  'VIC':'Victoria','Victoria':'Victoria',
  'TAS':'Tasmania','Tasmania':'Tasmania',
  'SA':'South Australia','South Australia':'South Australia',
  'WA':'Western Australia','Western Australia':'Western Australia',
  'NT':'Northern Territory','Northern Territory':'Northern Territory',
  'ACT':'Australian Capital Territory','Australian Capital Territory':'Australian Capital Territory'
};

var NDVI_BASELINE = 2005;
var NDVI_COMPARISON = 2024;
var DELTA_THR = 0.05; // absolute |ΔNDVI| threshold

/*  DATA  */
var mangroves   = ee.FeatureCollection(ASSETS.mangroves);
var seagrasses  = ee.FeatureCollection(ASSETS.seagrasses);
var saltmarshes = ee.FeatureCollection(ASSETS.saltmarshes);
var aoi = mangroves.merge(seagrasses).merge(saltmarshes);

var ausStates = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME','Australia'));
var ausGeom = ausStates.geometry();
var ausMask = ee.Image().byte().paint({featureCollection: ausStates, color:1}).selfMask();

// Ecosystem masks (byte; self-masked)
var mMask  = ee.Image().byte().paint({featureCollection: mangroves,   color:1}).selfMask();
var sgMask = ee.Image().byte().paint({featureCollection: seagrasses,  color:1}).selfMask();
var smMask = ee.Image().byte().paint({featureCollection: saltmarshes, color:1}).selfMask();
var bcAllMask = mMask.unmask().add(sgMask.unmask()).add(smMask.unmask()).gt(0).selfMask();

/*  NDVI (with light client-side caching)  */
var _ndviCache = {};
function ndviImageForYear(y){
  var key = String(y);
  if (_ndviCache[key]) return _ndviCache[key];
  var yNum = ee.Number(y);
  var img = ee.ImageCollection('MODIS/061/MOD13Q1')
              .filterDate(ee.Date.fromYMD(yNum,1,1), ee.Date.fromYMD(yNum.add(1),1,1))
              .select('NDVI')
              .median()
              .multiply(0.0001)
              .rename('NDVI')
              .updateMask(ausMask);
  _ndviCache[key] = img;
  return img;
}

function renderNdviForYear(y){
  var yNum = parseInt(y, 10); if (isNaN(yNum)) yNum = 2024;
  var ndvi = ndviImageForYear(yNum);
  lyrNDVI.setEeObject(ndvi);
  lyrNDVI.setVisParams({min:0, max:0.9, palette:PALETTES.ndvi});
  lyrNDVI.setName('NDVI (MODIS, ' + yNum + ')');
}

function renderDeltaBC(t1, t2, thrAbs){
  var diff = ndviImageForYear(t2).subtract(ndviImageForYear(t1)).updateMask(bcAllMask);
  if (thrAbs && thrAbs > 0) diff = diff.updateMask(diff.abs().gte(thrAbs));
  lyrNDVI.setEeObject(diff);
  lyrNDVI.setVisParams({min:-0.25, max:0.25, palette:PALETTES.dndvi});
  lyrNDVI.setName('ΔNDVI in Blue-Carbon ('+t2+' − '+t1+')');
}

/*  PERFORMANCE BASE  */
var CRS_METERS = 'EPSG:3577';     // GDA94 / Australian Albers
var SCALE_M = 250;                // Match MODIS 250 m

var pxArea_m2 = ee.Image.pixelArea().reproject({crs: CRS_METERS, scale: SCALE_M});
var haImgFix  = pxArea_m2.divide(10000);

var mMaskFix   = mMask.reproject({crs: CRS_METERS, scale: SCALE_M});
var sgMaskFix  = sgMask.reproject({crs: CRS_METERS, scale: SCALE_M});
var smMaskFix  = smMask.reproject({crs: CRS_METERS, scale: SCALE_M});
var allMaskFix = bcAllMask.reproject({crs: CRS_METERS, scale: SCALE_M});

function deltaNdvi(y1, y2){
  return ndviImageForYear(y2)
           .subtract(ndviImageForYear(y1))
           .rename('dNDVI')
           .updateMask(ausMask)
           .reproject({crs: CRS_METERS, scale: SCALE_M});
}

/*  FAST STATISTICS  */
function areaHaBands() {
  return ee.Image.cat([
    haImgFix.updateMask(mMaskFix).rename('ha_m'),
    haImgFix.updateMask(sgMaskFix).rename('ha_sg'),
    haImgFix.updateMask(smMaskFix).rename('ha_sm'),
    haImgFix.updateMask(allMaskFix).rename('ha_all')
  ]);
}

function fastStatsDict(geom, baseline, comparison, thrAbs) {
  var d = deltaNdvi(baseline, comparison);

  var lossMaskFix = d.lte(ee.Number(thrAbs).multiply(-1));
  var gainMaskFix = d.gte(ee.Number(thrAbs));

  var haBands = areaHaBands();
  var lossHa  = haBands.updateMask(lossMaskFix);
  var gainHa  = haBands.updateMask(gainMaskFix);

  var stack = ee.Image.cat([
    haBands.rename(['ha_m','ha_sg','ha_sm','ha_all']),
    lossHa .rename(['loss_m','loss_sg','loss_sm','loss_all']),
    gainHa .rename(['gain_m','gain_sg','gain_sm','gain_all'])
  ]);

  var gEff = ee.Geometry(geom).intersection(allMaskFix.geometry(), 1);

  var sums = stack.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: gEff,
    scale: SCALE_M,
    maxPixels: 1e13,
    tileScale: 4,
    bestEffort: true
  });

  return ee.Dictionary({
    areas: ee.Dictionary({
      m:  sums.get('ha_m'),
      sg: sums.get('ha_sg'),
      sm: sums.get('ha_sm'),
      all:sums.get('ha_all')
    }),
    loss: ee.Dictionary({
      m:  sums.get('loss_m'),
      sg: sums.get('loss_sg'),
      sm: sums.get('loss_sm'),
      all:sums.get('loss_all')
    }),
    gain: ee.Dictionary({
      m:  sums.get('gain_m'),
      sg: sums.get('gain_sg'),
      sm:  sums.get('gain_sm'),
      all: sums.get('gain_all')
    })
  });
}

/*  GEOMETRY HELPERS & OUTLINES  */
function toGaulName(name){ return GAUL_NAME_MAP[name] || name; }
function geomForState(stateVal){
  return (stateVal === 'All States')
    ? ausGeom
    : ee.Feature(ausStates.filter(ee.Filter.eq('ADM1_NAME', toGaulName(stateVal))).first()).geometry();
}

var OUTLINE_COLOR = '#FFFFFF';
function outlineImage(fc, width){
  return ee.Image().byte().paint({featureCollection: fc, color:1, width:width})
          .visualize({palette:[OUTLINE_COLOR], opacity:1});
}
var allStatesOutlineImg = outlineImage(ausStates, 1);

/*  MAP SETUP  */
var map = ui.Map();
map.setOptions('HYBRID');
map.setCenter(134.5, -25.0, 4);
map.setControlVisibility({layerList:true, mapTypeControl:true});

// Layers (preserve order)
var lyrSelectedOutline = ui.Map.Layer(ee.Image(0), {}, 'Selected State', false);
var lyrAllStatesOutline = ui.Map.Layer(allStatesOutlineImg, {}, 'All States (reference)', true);

var lyrMangrove  = ui.Map.Layer(mangroves.style(STYLES.mangrove),   {}, 'Mangroves',   false);
var lyrSeagrass  = ui.Map.Layer(seagrasses.style(STYLES.seagrass),  {}, 'Seagrasses',  false);
var lyrSaltmarsh = ui.Map.Layer(saltmarshes.style(STYLES.saltmarsh),{}, 'Saltmarshes', false);

var lyrNDVI = ui.Map.Layer(
  ee.Image(0), {min:0, max:0.9, palette:PALETTES.ndvi}, 'NDVI (MODIS, year)', false
);

// Insert in z-order: NDVI lowest, then habitats, then outlines
map.layers().reset([lyrNDVI, lyrMangrove, lyrSeagrass, lyrSaltmarsh, lyrAllStatesOutline, lyrSelectedOutline]);

/* ========== UI HELPERS ========== */
function hr(){ return ui.Panel([], null, {height:'1px', backgroundColor:THEME.border, margin:'8px 0'}); }
function heading(txt){ return ui.Label(txt, {fontWeight:'bold', color:THEME.text, margin:'8px 0 6px 0'}); }
function small(txt){ return ui.Label(txt, {color:THEME.subText, margin:'0 0 6px 0'}); }

function makeCard(title){
  var card = ui.Panel({style:{backgroundColor:THEME.cardBg, border:'1px solid '+THEME.border, margin:'8px', padding:'0'}});
  var header = ui.Panel([ ui.Label(title, {fontWeight:'bold', fontSize:'15px', color:THEME.text, margin:'0', padding:'8px 10px'}) ],
                        ui.Panel.Layout.flow('horizontal'),
                        {backgroundColor:THEME.cardHdr, padding:'0', margin:'0'});
  var body = ui.Panel({style:{padding:'8px'}});
  card.add(header); card.add(body);
  return {card:card, body:body};
}

/*  TIME SERIES (Inspector)  */
var yearsFull = ee.List.sequence(2001, 2024);
function ndviCollectionAnnual(){
  return ee.ImageCollection.fromImages(
    yearsFull.map(function(y){
      y = ee.Number(y);
      return ndviImageForYear(y)
        .set('system:time_start', ee.Date.fromYMD(y,7,1).millis())
        .set('year', y);
    })
  );
}
function ndviTimeSeriesAtPoint(pt){
  var col = ndviCollectionAnnual().select('NDVI');
  return ui.Chart.image.series({
      imageCollection: col, region: pt, reducer: ee.Reducer.mean(),
      scale: 250, xProperty: 'system:time_start'
    })
    .setChartType('LineChart')
    .setOptions({
      title: 'NDVI time series (250 m)',
      hAxis: {title: 'Year'}, vAxis: {title: 'NDVI'},
      legend: {position:'none'}, lineWidth:2, pointSize:3
    });
}

/*  SIDEBAR UI  */
var sidebar = ui.Panel({layout: ui.Panel.Layout.flow('vertical'),
  style:{width:'360px', padding:'0', backgroundColor:'rgba(0,0,0,0)', margin:'0'}});

var cardIntro = makeCard('Australia Blue-Carbon Explorer');
cardIntro.body.add(small('Toggle ecosystems, filter by state, and view compact statistics.'));

var cardLayers = makeCard('Layers');
var cbMangrove  = ui.Checkbox({label:'Show Mangroves',  value:false, onChange:function(v){ lyrMangrove.setShown(v); }});
var cbSeagrass  = ui.Checkbox({label:'Show Seagrasses', value:false, onChange:function(v){ lyrSeagrass.setShown(v); }});
var cbSaltmarsh = ui.Checkbox({label:'Show Saltmarshes',value:false, onChange:function(v){ lyrSaltmarsh.setShown(v); }});
cardLayers.body.add(cbMangrove); cardLayers.body.add(cbSeagrass); cardLayers.body.add(cbSaltmarsh);

var cardState = makeCard('Filter by state');
var stateContainer = ui.Panel({style:{margin:'4px 0 0 0'}});
stateContainer.add(small('Detecting fields and loading states…'));
cardState.body.add(stateContainer);

var cardInspector = makeCard('Inspector (click map)');
var chartHolder = ui.Panel();
cardInspector.body.add(chartHolder);

var cardStats = makeCard('State-wise statistics');
var statsPanel = ui.Panel();
cardStats.body.add(statsPanel);

sidebar.widgets().reset([cardIntro.card, cardLayers.card, cardState.card, cardInspector.card, cardStats.card]);

/* ========== RIGHT FLOAT UI ========== */
var quick = ui.Panel({
  style:{position:'top-right', padding:'10px', margin:'8px', backgroundColor:THEME.cardBg, border:'1px solid '+THEME.border}
});
quick.add(ui.Label('Map Controls', {fontWeight:'bold', color:THEME.text, margin:'0 0 6px 0'}));
quick.add(ui.Checkbox({label:'All-states outline', value:true, onChange:function(v){ lyrAllStatesOutline.setShown(v); }}));

quick.add(heading('Raster'));
var ndviYear = ui.Select({items:YEARS, value:'2024', style:{stretch:'horizontal'}});
var ndviToggle = ui.Checkbox({label:'NDVI (MODIS, annual)', value:false, onChange:function(v){ lyrNDVI.setShown(v); }});
var ndviAlpha = ui.Slider({min:0, max:1, value:0.35, step:0.05, onChange:function(v){ lyrNDVI.setOpacity(v); }});
quick.add(ndviYear);
quick.add(ndviToggle);
quick.add(ui.Label('NDVI opacity', {margin:'6px 0 2px 0', color:THEME.subText}));
quick.add(ndviAlpha);

/* Legend */
var legend = ui.Panel({
  style:{position:'bottom-left', padding:'8px', margin:'12px', backgroundColor:THEME.cardBg, border:'1px solid '+THEME.border}
});
legend.add(ui.Label('Legend', {fontWeight:'bold', color:THEME.text}));
legend.add(ui.Label('■ Mangroves',  {color:'#1b9e77', margin:'2px 0'}));
legend.add(ui.Label('■ Seagrasses', {color:'#377eb8', margin:'2px 0'}));
legend.add(ui.Label('■ Saltmarshes',{color:'#e6ab02', margin:'2px 0'}));

/*  Categorical NDVI legend  */
var NDVI_CLASSES = {
  breaks: [0.00, 0.20, 0.45, 0.65, 0.90],
  colors: ['#d73027', '#fdae61', '#b8e186', '#1f78b4'],
  labels: ['Bare / Water', 'Low vegetation', 'Moderate vegetation', 'Dense vegetation']
};
var _NDVI_CAT_WIDGETS = [];
function _ndviCatLegendClear(){ _NDVI_CAT_WIDGETS.forEach(function(w){ legend.remove(w); }); _NDVI_CAT_WIDGETS = []; }
function _ndviCatLegendAdd(w){ _NDVI_CAT_WIDGETS.push(w); legend.add(w); }
function setNdviCategoricalLegend() {
  _ndviCatLegendClear();
  _ndviCatLegendAdd(ui.Label('NDVI (current year)', {fontWeight:'bold', color: THEME.text, margin:'8px 0 6px 0'}));
  function row(hex, text){
    var p = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style:{margin:'0 0 4px 0'}});
    p.add(ui.Label('■', {color: hex, margin:'0 6px 0 0', fontWeight:'bold'}));
    p.add(ui.Label(text, {color: THEME.text}));
    _ndviCatLegendAdd(p);
  }
  row(NDVI_CLASSES.colors[0], NDVI_CLASSES.labels[0]);
  row(NDVI_CLASSES.colors[1], NDVI_CLASSES.labels[1]);
  row(NDVI_CLASSES.colors[2], NDVI_CLASSES.labels[2]);
  row(NDVI_CLASSES.colors[3], NDVI_CLASSES.labels[3]);
}
/*  categorical legend  */

/* Layout */
ui.root.clear();
ui.root.add(ui.SplitPanel({firstPanel: sidebar, secondPanel: map, orientation:'horizontal', wipe:false, style:{stretch:'both'}}));
map.add(legend);
map.add(quick);

/* Show the categorical NDVI legend immediately */
setNdviCategoricalLegend();

/*  INTERACTION LOGIC  */
var stateSelect, quickState, applyState;
function refreshYear(y){
  renderNdviForYear(y);
  if (applyState && stateSelect) applyState(stateSelect.getValue());
}
refreshYear('2024');

/* Field detection */
function detectField(fc, candidates, cb){
  fc.first().propertyNames().evaluate(function(names){
    names = names || [];
    var found = null;
    for (var i=0; i<candidates.length; i++){
      if (names.indexOf(candidates[i]) !== -1){ found = candidates[i]; break; }
    }
    cb(found);
  });
}
function distinctValues(fc, field, cb){
  if (!field){ cb([]); return; }
  fc.aggregate_array(field).distinct().sort().evaluate(function(vals){
    vals = (vals || []).filter(function(x){
      return x !== null && x !== undefined && x !== '' && x !== 'null' && x !== 'None';
    });
    cb(vals);
  });
}

/* Build state dropdowns after detecting fields & distinct values */
var mStateField, sgStateField, smStateField;
detectField(mangroves, STATE_CANDIDATES, function(msf){
  mStateField = msf;
  detectField(seagrasses, STATE_CANDIDATES, function(ssf){
    sgStateField = ssf;
    detectField(saltmarshes, STATE_CANDIDATES, function(slf){
      smStateField = slf;

      distinctValues(mangroves, mStateField, function(mStates){
        distinctValues(seagrasses, sgStateField, function(sgStates){
          distinctValues(saltmarshes, smStateField, function(smStates){

            var unionSet = {};
            (mStates||[]).forEach(function(x){ unionSet[x]=true; });
            (sgStates||[]).forEach(function(x){ unionSet[x]=true; });
            (smStates||[]).forEach(function(x){ unionSet[x]=true; });

            var items = ['All States'].concat(Object.keys(unionSet).sort());

            stateContainer.clear();
            stateSelect = ui.Select({items:items, value:'All States', placeholder:'Filter by state', style:{width:'100%'}});
            stateContainer.add(stateSelect);

            quick.add(ui.Label('State', {margin:'8px 0 2px 0', color:THEME.text, fontWeight:'bold'}));
            quickState = ui.Select({items:items, value:'All States', style:{stretch:'horizontal'}});
            quick.add(quickState);

            /* Main state apply */
            applyState = function(stateVal){
              var isAll = (stateVal === 'All States');

              var m  = (!mStateField || isAll) ? mangroves   : mangroves .filter(ee.Filter.eq(mStateField,  stateVal));
              var sg = (!sgStateField|| isAll) ? seagrasses  : seagrasses.filter(ee.Filter.eq(sgStateField, stateVal));
              var sm = (!smStateField|| isAll) ? saltmarshes : saltmarshes.filter(ee.Filter.eq(smStateField, stateVal));

              lyrMangrove.setEeObject(m.style(STYLES.mangrove));
              lyrSeagrass.setEeObject(sg.style(STYLES.seagrass));
              lyrSaltmarsh.setEeObject(sm.style(STYLES.saltmarsh));

              // Selected outline + view
              if (!isAll){
                var sFC = ausStates.filter(ee.Filter.eq('ADM1_NAME', toGaulName(stateVal)));
                lyrSelectedOutline.setEeObject(outlineImage(sFC, 3));
                lyrSelectedOutline.setShown(true);
                map.centerObject(ee.Feature(sFC.first()).geometry().centroid(1000), 7);
              } else {
                lyrSelectedOutline.setShown(false);
                map.setCenter(134.5, -25.0, 4);
              }

              //  Stats
              statsPanel.clear();
              statsPanel.add(ui.Label('Computing statistics…', {color:THEME.subText, margin:'0 0 6px 0'}));

              var gSel = geomForState(stateVal);
              var baseY = NDVI_BASELINE;
              var compY = NDVI_COMPARISON;
              var thr   = DELTA_THR;

              fastStatsDict(gSel, baseY, compY, thr).evaluate(function(d){
                statsPanel.clear();
                if (!d){ statsPanel.add(ui.Label('No statistics returned. Try again.', {color:'#b00'})); return; }

                var A = d.areas || {};
                var L = d.loss  || {};
                var G = d.gain  || {};

                function fmt(n,dg){ return (n==null ? '—' : Number(n).toLocaleString(undefined,{maximumFractionDigits:dg||1})); }
                function row(label,c2,c3,bold){
                  return ui.Panel([
                    ui.Label(label,{width:'140px', margin:'0 8px 0 0', color:THEME.text, fontWeight:bold?'bold':'normal'}),
                    ui.Label(c2,   {width:'120px', margin:'0 8px 0 0', color:THEME.text, fontWeight:bold?'bold':'normal'}),
                    ui.Label(c3,   {color:THEME.text, fontWeight:bold?'bold':'normal'})
                  ], ui.Panel.Layout.flow('horizontal'));
                }

                // Area table (km² + ha)
                statsPanel.add(heading('Habitat area'));
                var head = ui.Panel({style:{backgroundColor:THEME.cardHdr, padding:'6px 6px 6px 0', margin:'0 0 6px 0'}});
                head.add(row('Ecosystem','Area (km²)','Area (ha)',true));
                statsPanel.add(head);

                var toKm2 = function(ha){ return (ha==null) ? null : Number(ha)/100; };
                statsPanel.add(row('Mangroves',   fmt(toKm2(A.m)),   fmt(A.m)));
                statsPanel.add(row('Seagrasses',  fmt(toKm2(A.sg)),  fmt(A.sg)));
                statsPanel.add(row('Saltmarshes', fmt(toKm2(A.sm)),  fmt(A.sm)));
                statsPanel.add(hr());
                statsPanel.add(row('All habitats', fmt(toKm2(A.all)), fmt(A.all), true));

                // ΔNDVI Loss/Gain as % area
                function pct(valueHa, denomHa){
                  if (valueHa==null || denomHa==null || Number(denomHa)===0) return '—';
                  return (100*Number(valueHa)/Number(denomHa)).toFixed(1)+'%';
                }
                function pctLine(label, lossP, gainP){
                  var line = ui.Panel([
                    ui.Label(label+' — ', {color:THEME.text, margin:'0 8px 0 0'}),
                    ui.Label('Loss: '+lossP, {color:'#b91c1c', margin:'0 10px 0 0'}),
                    ui.Label('Gain: '+gainP, {color:'#15803d'})
                  ], ui.Panel.Layout.flow('horizontal'));
                  statsPanel.add(line);
                }

                statsPanel.add(hr());
                statsPanel.add(heading('ΔNDVI (Normalised Difference Vegetation Index, NDVI) Loss/Gain (% |Δ| ≥ '+thr+') — '+stateVal+' ['+compY+' − '+baseY+']'));

                pctLine('All habitats', pct(L.all, A.all), pct(G.all, A.all));
                pctLine('Mangroves',    pct(L.m,   A.m),   pct(G.m,   A.m));
                pctLine('Seagrasses',   pct(L.sg,  A.sg),  pct(G.sg,  A.sg));
                pctLine('Saltmarshes',  pct(L.sm,  A.sm),  pct(G.sm,  A.sm));

                // CSV export (areas)
                statsPanel.add(hr());
                statsPanel.add(ui.Button({
                  label:'Export areas CSV ('+stateVal+')',
                  onClick:function(){
                    var fc = ee.FeatureCollection([
                      ee.Feature(null, {state:stateVal, ecosystem:'Mangroves',   area_km2:toKm2(A.m),   area_ha:A.m}),
                      ee.Feature(null, {state:stateVal, ecosystem:'Seagrasses',  area_km2:toKm2(A.sg),  area_ha:A.sg}),
                      ee.Feature(null, {state:stateVal, ecosystem:'Saltmarshes', area_km2:toKm2(A.sm),  area_ha:A.sm}),
                      ee.Feature(null, {state:stateVal, ecosystem:'All habitats',area_km2:toKm2(A.all), area_ha:A.all})
                    ]);
                    var fname = 'bluecarbon_area_by_ecosystem_' + stateVal.replace(/\s+/g,'_');
                    Export.table.toDrive({collection:fc, description:fname, fileNamePrefix:fname, fileFormat:'CSV'});
                    print('Created an export task for:', fname);
                  },
                  style:{margin:'6px 0'}
                }));
              });
            };

            // Sync dropdowns
            stateSelect.onChange(function(v){ quickState.setValue(v, false); applyState(v); });
            quickState.onChange(function(v){ stateSelect.setValue(v, false); applyState(v); });

            // Year change -> refresh raster + stats
            ndviYear.onChange(function(y){ refreshYear(y); });

            // Initial draw
            applyState('All States');

            print('Detected state fields:', {mangroves:{state:mStateField}, seagrasses:{state:sgStateField}, saltmarshes:{state:smStateField}});
          });
        });
      });
    });
  });
});

/*  MAP CLICK (Inspector only: NDVI time series)  */
map.onClick(function(coords){
  chartHolder.clear();
  var pt = ee.Geometry.Point([coords.lon, coords.lat]);
  chartHolder.add(ndviTimeSeriesAtPoint(pt));
});

