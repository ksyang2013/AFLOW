// aflowlib_webapp_bands.cpp automatic generated
std::string AFLOW_WEBAPP_BANDS_JS="\
// aflowlib_webapp_bands.js \n\
// \n\
// author: Harvey Shi \n\
// edited: Geena Gomez (February 2018) \n\
// edited: Pauline Colinet (May 2018) \n\
// \n\
// JS file for the interactive bands plot and \n\
//    density of states plot on the entry page. \n\
 \n\
/* \n\
function getMax(arrayOfArray) { \n\
    for (var i = 1; i< arrayOfArray.length; i++) { \n\
      arr = arrayOfArray[i]; \n\
      intMax = 0; \n\
      var max = arr.reduce(function(a, b) { \n\
            return Math.max(a, b); \n\
      }); \n\
      finalMax = Math.max(max, intMax); \n\
      intMax = max; \n\
      return finalMax \n\
    } \n\
}; \n\
*/ \n\
 \n\
$(document).ready(function() { \n\
  var json = d3_bands_data, \n\
    lineLabels = json.kpoint_labels_html, \n\
    lineLocs = json.kpoint_positions, \n\
    title = json.title, \n\
    yExtrema = [-5, 5]; //json.Emin, json.Emax]; \n\
 \n\
  // Reformat line labels \n\
  for (var i = 0; i < lineLabels.length; ++i) { \n\
    var label = lineLabels[i]; \n\
    label = label.replace(\"G\", \"G\"); \n\
    label = label.replace(\"Sigm\", \"Î£\"); \n\
    label = label.replace(\"_1\", \"<tspan class='subScript' dy='3'>1</tspan>\"); \n\
    lineLabels[i] = label; \n\
  } \n\
 \n\
  var spinP = json.pDOS_data.spin_polarized; \n\
  if (spinP == true) { \n\
    // Spin Polarized Case \n\
    var bandData = [json.bands_data_majority, json.bands_data_minority], \n\
      xValues = function(d) { \n\
        return d[0]; \n\
      }, \n\
      bandX = [d3.min(bandData[0], xValues), d3.max(bandData[0], xValues)], \n\
      numBands = json.n_bands, \n\
      rows = numBands * 2, \n\
      cols = bandData[0].length, \n\
      maxBandPts = cols, \n\
      bandArray = new Array(); \n\
    for (var i = 0; i < rows; ++i) { \n\
      bandArray.push([]); \n\
    } \n\
    for (var i = 0; i < rows; ++i) { \n\
      for (var j = 0; j < cols; ++j) { \n\
        if (i < numBands) { \n\
          bandArray[i].push({ \n\
            x: bandData[0][j][0], \n\
            y: bandData[0][j][i + 1] \n\
          }); \n\
        } else { \n\
          bandArray[i].push({ \n\
            x: bandData[1][j][0], \n\
            y: bandData[1][j][i - numBands + 1] \n\
          }); \n\
        } \n\
      } \n\
    } \n\
    var bandXArray = new Array(); \n\
    for (var j = 0; j < cols; ++j) { \n\
      bandXArray.push(bandData[0][j][0]); \n\
    } \n\
    var dosData = [ \n\
        json.pDOS_data.energy, \n\
        json.pDOS_data.sum_s_majority, \n\
        json.pDOS_data.sum_p_majority, \n\
        json.pDOS_data.sum_d_majority, \n\
        json.pDOS_data.sum_s_minority, \n\
        json.pDOS_data.sum_p_minority, \n\
        json.pDOS_data.sum_d_minority \n\
      ], \n\
      dosNames = [ \n\
        \"s majority\", \n\
        \"p majority\", \n\
        \"d majority\", \n\
        \"s minority\", \n\
        \"p minority\", \n\
        \"d minority\" \n\
      ], \n\
      numDos = dosNames.length; \n\
    var dosArray = new Array(), \n\
      rows = dosData.length - 1, \n\
      cols = dosData[0].length; \n\
 \n\
    //PC -- 2018/08/09 -- Calculates x range for DOS plot \n\
    var minDosX = 1000; \n\
    var maxDosX = -1000; \n\
    for (var i = 0; i < dosData[0].length; i++) { \n\
      if ((dosData[0][i] > -5) && (dosData[0][i] < 5)) { \n\
      var arrTemp = []  \n\
      for (var j = 1; j < dosData.length; j++) { \n\
          arrTemp.push(dosData[j][i]); \n\
       } \n\
      var minArrTemp = Math.min.apply(null, arrTemp); \n\
      minDosX = Math.min(minDosX, minArrTemp); \n\
      var maxArrTemp = Math.max.apply(null, arrTemp); \n\
      maxDosX = Math.max(maxDosX, maxArrTemp); \n\
      } \n\
    } \n\
    maxDosX = maxDosX +1; \n\
    minDosX = minDosX -1; \n\
    finalMaxDosX = Math.max(Math.abs(minDosX),maxDosX); \n\
    for (var i = 0; i < rows; ++i) { \n\
      dosArray.push([]); \n\
    } \n\
    for (var i = 0; i < rows; ++i) { \n\
      for (var j = 0; j < cols; ++j) { \n\
        dosArray[i].push({ \n\
          x: dosData[i + 1][j], \n\
          y: dosData[0][j] \n\
        }); \n\
      } \n\
    } \n\
 \n\
    $(\"#bandLegend\").fadeIn(200); \n\
    $(\"#dosText\").css(\"height\", \"60px\"); \n\
    var dosYArray = dosData[0], \n\
      dosColor = [ \n\
        \"#66ff29\", \n\
        \"#0662F9\", \n\
        \"#FA7943\", \n\
        \"#66ff29\", \n\
        \"#0662F9\", \n\
        \"#FA7943\" \n\
      ], \n\
      legendNames = [\"s\", \"p\", \"d\"], \n\
      dosData = [ \n\
        json.pDOS_data.sum_s_majority, \n\
        json.pDOS_data.sum_p_majority, \n\
        json.pDOS_data.sum_d_majority, \n\
        json.pDOS_data.sum_s_minority, \n\
        json.pDOS_data.sum_p_minority, \n\
        json.pDOS_data.sum_d_minority \n\
      ]; \n\
    if (numDos == 4) { \n\
      dosData.push(json.pDOS_data.sum_f); \n\
    } \n\
 \n\
    var min = d3.min( \n\
      dosData.map(function(array) { \n\
        return d3.min(array); \n\
      }) \n\
    ); \n\
    var max = d3.max( \n\
      dosData.map(function(array) { \n\
        return d3.max(array); \n\
      }) \n\
    ); \n\
    // PC --  2018/08/09 \n\
  // dosX = [-10, 10]; \n\
    //dosX = [min-1, max+1]; \n\
     dosX = [-finalMaxDosX, finalMaxDosX]; \n\
 \n\
    // Transpose matrix \n\
    dosData = dosData[0].map(function(col, i) { \n\
      return dosData.map(function(row) { \n\
        return row[i]; \n\
      }); \n\
    }); \n\
    if (numDos == 4) { \n\
      dosData.push(json.pDOS_data.sum_f); \n\
    } \n\
  } else { \n\
    // Not Spin Polarized Case \n\
    var bandData = json.bands_data, \n\
      numBands = json.n_bands, \n\
      xValues = function(d) { \n\
        return d[0]; \n\
      }, \n\
      bandX = [d3.min(bandData, xValues), d3.max(bandData, xValues)], \n\
      bandArray = new Array(), \n\
      rows = numBands, \n\
      cols = bandData.length, \n\
      maxBandPts = cols; \n\
    for (var i = 0; i < rows; ++i) { \n\
      bandArray.push([]); \n\
    } \n\
    for (var i = 0; i < rows; ++i) { \n\
      for (var j = 0; j < cols; ++j) { \n\
        bandArray[i].push({ \n\
          x: bandData[j][0], \n\
          y: bandData[j][i + 1] \n\
        }); \n\
      } \n\
    } \n\
    var bandXArray = new Array(); \n\
    for (var j = 0; j < cols; ++j) { \n\
      bandXArray.push(bandData[j][0]); \n\
    } \n\
    var dosData = [ \n\
        json.pDOS_data.energy, \n\
        json.pDOS_data.sum_s, \n\
        json.pDOS_data.sum_p, \n\
        json.pDOS_data.sum_d \n\
      ], \n\
      dosNames = json.pDOS_data.orbitals, \n\
      numDos = dosNames.length; \n\
    if (numDos == 4) { \n\
      dosData.push(json.pDOS_data.sum_f); \n\
    } \n\
    var dosArray = new Array(), \n\
      rows = dosData.length - 1, \n\
      cols = dosData[0].length; \n\
    //PC -- 2018/08/09 -- Calculates x range for DOS plot \n\
    var maxDosX = -1000; \n\
    for (var i = 0; i < dosData[0].length; i++) { \n\
      if ((dosData[0][i] > -5) && (dosData[0][i] < 5)) { \n\
      var arrTemp = []  \n\
      for (var j = 1; j < dosData.length; j++) { \n\
          arrTemp.push(dosData[j][i]); \n\
       } \n\
      var maxArrTemp = Math.max.apply(null, arrTemp); \n\
      maxDosX = Math.max(maxDosX, maxArrTemp); \n\
      } \n\
    } \n\
    maxDosX = maxDosX +1; \n\
 \n\
    for (var i = 0; i < rows; ++i) { \n\
      dosArray.push([]); \n\
    } \n\
    for (var i = 0; i < rows; ++i) { \n\
      for (var j = 0; j < cols; ++j) { \n\
        dosArray[i].push({ \n\
          x: dosData[i + 1][j], \n\
          y: dosData[0][j] \n\
        }); \n\
      } \n\
    } \n\
 \n\
    var dosYArray = dosData[0], \n\
      dosColor = [\"#66ff29\", \"#0662F9\", \"#FA7943\", \"#19e3e3\"], \n\
      dosData = [ \n\
        json.pDOS_data.sum_s, \n\
        json.pDOS_data.sum_p, \n\
        json.pDOS_data.sum_d \n\
      ]; \n\
    if (numDos == 4) { \n\
      dosData.push(json.pDOS_data.sum_f); \n\
      var dosData3 = dosData[3], \n\
        dosDataFiltered3 = new Array(); \n\
    } \n\
 \n\
    var dosData0 = dosData[0], \n\
      dosData1 = dosData[1], \n\
      dosData2 = dosData[2]; // separating dosData data sets \n\
 \n\
    var rangeIndexArray = new Array(); \n\
    for (i = 0; i < cols; ++i) { \n\
      if (dosYArray[i] <= 10 && dosYArray[i] >= -10) { \n\
        rangeIndexArray.push(i); \n\
      } \n\
    } \n\
    var rangeStart = rangeIndexArray[0], \n\
      rangeEnd = rangeIndexArray.length + rangeStart; // beginning and end of the range of relevant energy value indexes \n\
 \n\
    var dosDataFiltered0 = new Array(), \n\
      dosDataFiltered1 = new Array(), \n\
      dosDataFiltered2 = new Array(); \n\
    for (i = rangeStart; i < rangeEnd; ++i) { \n\
      dosDataFiltered0.push(dosData0[i]); \n\
    } \n\
    for (i = rangeStart; i < rangeEnd; ++i) { \n\
      dosDataFiltered1.push(dosData1[i]); \n\
    } \n\
    for (i = rangeStart; i < rangeEnd; ++i) { \n\
      dosDataFiltered2.push(dosData2[i]); \n\
    } \n\
    if (numDos == 4) { \n\
      for (i = rangeStart; i < rangeEnd; ++i) { \n\
        dosDataFiltered3.push(dosData3[i]); \n\
      } \n\
      var dosDataFiltered = [ \n\
        dosDataFiltered0, \n\
        dosDataFiltered1, \n\
        dosDataFiltered2, \n\
        dosDataFiltered3 \n\
      ]; // new dosData array with values within relevant indexes \n\
    } else { \n\
      var dosDataFiltered = [ \n\
        dosDataFiltered0, \n\
        dosDataFiltered1, \n\
        dosDataFiltered2 \n\
      ]; \n\
    } \n\
 \n\
    //for (var i = 0; i< dosData.length; i++) {   \n\
    //} \n\
 \n\
    var max = d3.max( \n\
      dosDataFiltered.map(function(array) { \n\
        return d3.max(array); \n\
      }) \n\
    ); \n\
    // PC -- 2018/08/09  \n\
    //dosX = [0, max + 2]; \n\
     dosX = [0, maxDosX]; \n\
 \n\
    // [OBSOLETE]    var max = d3.max(dosData.map(function (array) { \n\
    // [OBSOLETE]        return d3.max(array); \n\
    // [OBSOLETE]   })); \n\
 \n\
    // Transpose matrix \n\
    dosData = dosData[0].map(function(col, i) { \n\
      return dosData.map(function(row) { \n\
        return row[i]; \n\
      }); \n\
    }); \n\
  } \n\
  var bands = d3.select(\"#bands_wrapper\"), \n\
    WIDTH = 900, \n\
    HEIGHT = 550, \n\
    MARGINS = { \n\
      top: 30, \n\
      right: 10, \n\
      bottom: 30, \n\
      left: 60 \n\
    }, \n\
    d3xR = d3.scale \n\
      .linear() \n\
      .range([MARGINS.left, WIDTH - MARGINS.right]) \n\
      .domain([bandX[0], bandX[1]]), \n\
    d3yR = d3.scale \n\
      .linear() \n\
      .range([HEIGHT - MARGINS.top, MARGINS.bottom]) \n\
      .domain([yExtrema[0], yExtrema[1]]), \n\
    xAxis = d3.svg \n\
      .axis() \n\
      .scale(d3xR) \n\
      .tickSubdivide(true) \n\
      .tickValues(false), \n\
    yAxis = d3.svg \n\
      .axis() \n\
      .scale(d3yR) \n\
      .tickSize(5) \n\
      .orient(\"left\") \n\
      .tickSubdivide(true); \n\
  var dosOffset = WIDTH; \n\
  bands \n\
    .append(\"clipPath\") \n\
    .attr(\"id\", \"clip\") \n\
    .append(\"rect\") \n\
    .attr(\"width\", WIDTH - MARGINS.left - MARGINS.right) \n\
    .attr(\"height\", HEIGHT - MARGINS.top - MARGINS.bottom) \n\
    .attr(\"x\", MARGINS.left) \n\
    .attr(\"y\", MARGINS.top); \n\
 \n\
  bands \n\
    .append(\"svg:g\") \n\
    .attr(\"class\", \"x axis\") \n\
    .attr(\"transform\", \"translate(0,\" + (HEIGHT - MARGINS.bottom) + \")\") \n\
    .call(xAxis); \n\
  // Add k-point labels \n\
  for (var i = 0; i < lineLabels.length; ++i) { \n\
    var xPos; \n\
    if (i == lineLabels.length - 1) { \n\
      xPos = d3xR(d3.max(bandData[0], xValues)); \n\
    } else { \n\
      xPos = d3xR(lineLocs[i]); \n\
    } \n\
    var group = bands \n\
      .append(\"g\") \n\
      .attr(\"class\", \"lineLabels\") \n\
      .attr(\"data-xPos\", xPos) \n\
      .attr(\"transform\", function() { \n\
        return \"translate(\" + xPos + \",0)\"; \n\
      }); \n\
 \n\
    group \n\
      .append(\"text\") \n\
      .html(lineLabels[i]) \n\
      .attr(\"x\", 0) \n\
      .attr(\"y\", HEIGHT - MARGINS.bottom + 20) \n\
      .attr(\"text-anchor\", \"middle\"); \n\
 \n\
    group \n\
      .append(\"svg:line\") \n\
      .attr(\"class\", \"line\") \n\
      .attr(\"stroke-width\", 1) \n\
      .attr(\"x1\", 0) \n\
      .attr(\"y1\", MARGINS.top) \n\
      .attr(\"x2\", 0) \n\
      .attr(\"y2\", HEIGHT - MARGINS.top); \n\
  } \n\
 \n\
  // Pseudo-clipping mask for k-point lines and labels \n\
  bands \n\
    .append(\"rect\") \n\
    .attr(\"fill\", \"white\") \n\
    .attr(\"width\", MARGINS.left) \n\
    .attr(\"height\", HEIGHT - MARGINS.bottom) \n\
    .attr(\"x\", 0) \n\
    .attr(\"y\", 0); \n\
  bands \n\
    .append(\"rect\") \n\
    .attr(\"fill\", \"white\") \n\
    .attr(\"width\", MARGINS.left - 5) \n\
    .attr(\"height\", MARGINS.bottom) \n\
    .attr(\"x\", 0) \n\
    .attr(\"y\", HEIGHT - MARGINS.bottom); \n\
  bands \n\
    .append(\"rect\") \n\
    .attr(\"fill\", \"white\") \n\
    .attr(\"width\", MARGINS.right) \n\
    .attr(\"height\", HEIGHT - MARGINS.bottom) \n\
    .attr(\"x\", WIDTH - MARGINS.right) \n\
    .attr(\"y\", 0); \n\
  bands \n\
    .append(\"rect\") \n\
    .attr(\"fill\", \"white\") \n\
    .attr(\"width\", MARGINS.right - 5) \n\
    .attr(\"height\", MARGINS.bottom) \n\
    .attr(\"x\", WIDTH - MARGINS.right + 5) \n\
    .attr(\"y\", HEIGHT - MARGINS.bottom); \n\
 \n\
  bands \n\
    .append(\"svg:g\") \n\
    .attr(\"class\", \"y axis\") \n\
    .attr(\"transform\", \"translate(\" + MARGINS.left + \",0)\") \n\
    .call(yAxis); \n\
 \n\
  bands \n\
    .append(\"text\") \n\
    .attr(\"class\", \"titles\") \n\
    .attr(\"x\", WIDTH / 2) \n\
    .attr(\"y\", 15) \n\
    .attr(\"text-anchor\", \"middle\") \n\
    .text(title); \n\
 \n\
  bands \n\
    .append(\"text\") \n\
    .attr(\"class\", \"titles\") \n\
    .attr(\"x\", -(HEIGHT / 2 + MARGINS.top + 20)) \n\
    .attr(\"y\", 0) \n\
    .attr(\"dy\", \".75em\") \n\
    .attr(\"transform\", \"rotate(-90)\") \n\
    .text(\"Energy (eV)\"); \n\
 \n\
  bands \n\
    .append(\"svg:svg\") \n\
    .attr(\"id\", \"bandMove\") \n\
    .attr(\"fill\", \"transparent\") \n\
    .attr(\"width\", WIDTH - MARGINS.left - MARGINS.right) \n\
    .attr(\"height\", HEIGHT - MARGINS.top - MARGINS.bottom) \n\
    .attr(\"x\", MARGINS.left) \n\
    .attr(\"y\", MARGINS.top); \n\
 \n\
  var line = d3.svg \n\
    .line() \n\
    .interpolate(\"linear\") \n\
    .x(function(d) { \n\
      return d3xR(d.x); \n\
    }) \n\
    .y(function(d) { \n\
      return d3yR(d.y); \n\
    }); \n\
 \n\
  if (spinP == true) { \n\
    var minColor = \"#ff0d00\"; \n\
    $(\".legendLine.min\").css(\"border-color\", \"#ff0d00\"); \n\
    var band_lines = bands \n\
      .selectAll(\".band_lines\") \n\
      .data(bandArray) \n\
      .enter() \n\
      .append(\"path\") \n\
      .attr(\"class\", \"band_lines\") \n\
      .attr(\"stroke\", function(d, i) { \n\
        if (i < numBands) { \n\
          return \"black\"; \n\
        } else { \n\
          return minColor; \n\
        } \n\
      }) \n\
      .attr(\"id\", function(d, i) { \n\
        if (i < numBands) { \n\
          return \"majSpin\"; \n\
        } else { \n\
          return \"minSpin\"; \n\
        } \n\
      }) \n\
      .attr(\"stroke-width\", 1) \n\
      .attr(\"fill\", \"none\") \n\
      .attr(\"clip-path\", \"url(#clip)\") \n\
      .attr(\"d\", line); \n\
 \n\
    var bandTotal = [], \n\
      max = bandData[0].length; \n\
    for (var i = 0; i < max; ++i) { \n\
      bandTotal.push(bandData[0][i].concat(bandData[1][i])); \n\
    } \n\
  } else { \n\
    var band_lines = bands.selectAll(\".band_lines\").data(bandArray); \n\
    band_lines \n\
      .enter() \n\
      .append(\"path\") \n\
      .attr(\"class\", \"band_lines\") \n\
      .attr(\"stroke\", \"rgba(0,0,0,.7)\") \n\
      .attr(\"stroke-width\", 1) \n\
      .attr(\"fill\", \"none\") \n\
      .attr(\"clip-path\", \"url(#clip)\") \n\
      .attr(\"d\", line); \n\
    var bandTotal = bandData; \n\
  } \n\
 \n\
  bands.on(\"mousemove\", mousemoveSpin).on(\"mouseover\", function() { \n\
    $(\"#dosTip\") \n\
      .attr(\"opacity\", 0) \n\
      .attr(\"cx\", 0 + \"px\") \n\
      .attr(\"cy\", 0 + \"px\"); \n\
    $(\"#dosText\").css(\"opacity\", 0); \n\
  }); \n\
  var bisect = d3.bisector(function(d) { \n\
    return d; \n\
  }).right; \n\
 \n\
  //modif Pauline \n\
  //Creation of two simple arrays (one for majority, one for minority) with the y coordinates of each band. \n\
  var bandMajority = []; \n\
  var bandMinority = []; \n\
  for (i = 0; i < bandTotal.length; i++) { \n\
    var bandM_i = []; \n\
    var bandm_i = []; \n\
    for (j = 1; j <= 0.5 * (bandTotal[0].length - 2); j++) { \n\
      bandM_i.push(bandTotal[i][j]); \n\
      bandm_i.push(bandTotal[i][j + 0.5 * bandTotal[0].length]); \n\
    } \n\
    bandMajority.push(bandM_i); \n\
    bandMinority.push(bandm_i); \n\
  } \n\
 \n\
  //modif Pauline: new mousemove function \n\
  //mousemove changed so that it takes into account the selected value of the Spin Selection button \n\
  //cf var spinSelectorValue \n\
  // \n\
  //modif Pauline: mousemoveSpin works also for non-spin polarized cases \n\
  //i.e. SpinP != true \n\
 \n\
  function mousemoveSpin() { \n\
    var mouseX = d3.mouse(this)[0], \n\
      mouseY = d3.mouse(this)[1]; \n\
    var xPos = d3xR.invert(mouseX), \n\
      xIndex = bisect(bandXArray, xPos); \n\
    if (xIndex == maxBandPts) { \n\
      --xIndex; \n\
    } \n\
 \n\
    var bandTot = []; \n\
    if (document.getElementById(\"spinBandsOptions\") == null) { \n\
      var spinSelectorValue = \"none\"; \n\
    } else { \n\
      var spinSelectorValue = document.getElementById(\"spinBandsOptions\").value; \n\
    } \n\
 \n\
    if (xPos > 0 && xPos < 0.01 + bandXArray[bandXArray.length - 1]) { \n\
      $(\"#bandTip\").attr(\"visibility\", \"visible\"); \n\
      if (spinSelectorValue == \"majority\") { \n\
        $(\"#bandTip rect\").attr(\"stroke\", \"black\"); \n\
        bandTot = bandMajority; \n\
        bandTotBis = bandMajority; \n\
      } else if (spinSelectorValue == \"minority\") { \n\
        bandTot = bandMinority; \n\
        bandTotBis = bandMinority; \n\
        $(\"#bandTip rect\").attr(\"stroke\", minColor); \n\
      } else if (spinSelectorValue == \"bothS\") { \n\
        bandTot = bandMajority; \n\
        bandTotBis = bandMinority; \n\
      } else if (spinSelectorValue == \"none\") { \n\
        bandTot = bandTotal; \n\
        bandTotBis = bandTotal; \n\
        $(\"#bandTip rect\").attr(\"stroke\", \"black\"); \n\
      } \n\
 \n\
      var yPos = d3yR.invert(mouseY), \n\
        i = bisect(bandTot[xIndex], yPos), \n\
        d0 = bandTot[xIndex][i - 1], \n\
        d1 = bandTot[xIndex][i], \n\
        index = yPos - d0 > d1 - yPos ? i : i - 1, \n\
        y1 = bandTot[Math.max(xIndex - 1, 0)][index], \n\
        y2 = bandTot[xIndex][index], \n\
        x1 = bandXArray[Math.max(xIndex - 1, 0)], \n\
        x2 = bandXArray[xIndex]; \n\
      if (x2 - x1 == 0) { \n\
        var yLocOne = 0; \n\
      } else { \n\
        var yLocOne = y1 + (xPos - x1) * (y2 - y1) / (x2 - x1); // Linear interpolation \n\
      } \n\
      var yPosBis = d3yR.invert(mouseY), \n\
        iBis = bisect(bandTotBis[xIndex], yPosBis), \n\
        d0Bis = bandTotBis[xIndex][iBis - 1], \n\
        d1Bis = bandTotBis[xIndex][iBis], \n\
        indexBis = yPosBis - d0Bis > d1Bis - yPosBis ? iBis : iBis - 1, \n\
        y1Bis = bandTotBis[Math.max(xIndex - 1, 0)][indexBis], \n\
        y2Bis = bandTotBis[xIndex][indexBis], \n\
        x1Bis = bandXArray[Math.max(xIndex - 1, 0)], \n\
        x2Bis = bandXArray[xIndex]; \n\
      if (x2Bis - x1Bis == 0) { \n\
        var yLocBis = 0; \n\
      } else { \n\
        var yLocBis = \n\
          y1Bis + (xPos - x1Bis) * (y2Bis - y1Bis) / (x2Bis - x1Bis); // Linear interpolation \n\
      } \n\
 \n\
      yLoc = \n\
        Math.abs(yPosBis - yLocBis) > Math.abs(yPos - yLocOne) \n\
          ? yLocOne \n\
          : yLocBis; \n\
 \n\
      if (yLocOne !== yLocBis) { \n\
        if (yLoc === yLocOne) { \n\
          $(\"#bandTip rect\").attr(\"stroke\", \"black\"); \n\
        } else { \n\
          $(\"#bandTip rect\").attr(\"stroke\", minColor); \n\
        } \n\
      } \n\
      if (d3xR(xPos) > 735) { \n\
        $(\"#bandTip\") \n\
          .attr(\"x\", d3xR(xPos) - 108 + \"px\") \n\
          .attr(\"y\", d3yR(yLoc) - 16 + \"px\"); \n\
        $(\"#bandTip circle\").attr(\"cx\", \"104\"); \n\
        $(\"#bandTip polygon\").attr(\"points\", \"96, 16, 89, 10, 89, 22 \"); \n\
        $(\"#bandTip rect\").attr(\"x\", 2); \n\
        $(\"#bandTip text\").attr(\"x\", 45); \n\
      } else { \n\
        $(\"#bandTip\") \n\
          .attr(\"x\", d3xR(xPos) - 6 + \"px\") \n\
          .attr(\"y\", d3yR(yLoc) - 16 + \"px\"); \n\
        $(\"#bandTip circle\").attr(\"cx\", \"6\"); \n\
        $(\"#bandTip polygon\").attr(\"points\", \"15, 16, 22, 10, 22, 22 \"); \n\
        $(\"#bandTip rect\").attr(\"x\", 21); \n\
        $(\"#bandTip text\").attr(\"x\", 65); \n\
      } \n\
      $(\"#bandText\").text(yLoc.toFixed(4) + \" eV\"); \n\
      checkBounds(mouseX, mouseY, d3yR(yLoc) + 15); \n\
    } else { \n\
      $(\"#bandTip\").attr(\"visibility\", \"hidden\"); \n\
    } \n\
  } \n\
 \n\
 function mousemove() { \n\
    var mouseX = d3.mouse(this)[0], \n\
      mouseY = d3.mouse(this)[1]; \n\
    var xPos = d3xR.invert(mouseX), \n\
      xIndex = bisect(bandXArray, xPos); \n\
    if (xIndex == maxBandPts) { \n\
      --xIndex; \n\
    } \n\
 \n\
    var yPos = d3yR.invert(mouseY), \n\
      i = bisect(bandTotal[xIndex], yPos), \n\
      d0 = bandTotal[xIndex][i - 1], \n\
      d1 = bandTotal[xIndex][i], \n\
      index = yPos - d0 > d1 - yPos ? i : i - 1, \n\
      y1 = bandTotal[Math.max(xIndex - 1, 0)][index], \n\
      y2 = bandTotal[xIndex][index], \n\
      x1 = bandXArray[Math.max(xIndex - 1, 0)], \n\
      x2 = bandXArray[xIndex]; \n\
    if (x2 - x1 == 0) { \n\
      var yLoc = 0; \n\
    } else { \n\
      var yLoc = y1 + (xPos - x1) * (y2 - y1) / (x2 - x1); // Linear interpolation \n\
    } \n\
    if (index > numBands) { \n\
      $(\"#bandTip rect\").attr(\"stroke\", minColor); \n\
    } else { \n\
      $(\"#bandTip rect\").attr(\"stroke\", \"black\"); \n\
    } \n\
    if (d3xR(xPos) > 735) { \n\
      $(\"#bandTip\") \n\
        .attr(\"x\", d3xR(xPos) - 108 + \"px\") \n\
        .attr(\"y\", d3yR(yLoc) - 16 + \"px\"); \n\
      $(\"#bandTip circle\").attr(\"cx\", \"104\"); \n\
      $(\"#bandTip polygon\").attr(\"points\", \"96, 16, 89, 10, 89, 22 \"); \n\
      $(\"#bandTip rect\").attr(\"x\", 2); \n\
      $(\"#bandTip text\").attr(\"x\", 45); \n\
    } else { \n\
      $(\"#bandTip\") \n\
        .attr(\"x\", d3xR(xPos) - 6 + \"px\") \n\
        .attr(\"y\", d3yR(yLoc) - 16 + \"px\"); \n\
      $(\"#bandTip circle\").attr(\"cx\", \"6\"); \n\
      $(\"#bandTip polygon\").attr(\"points\", \"15, 16, 22, 10, 22, 22 \"); \n\
      $(\"#bandTip rect\").attr(\"x\", 21); \n\
      $(\"#bandTip text\").attr(\"x\", 65); \n\
    } \n\
    $(\"#bandText\").text(yLoc.toFixed(4) + \" eV\"); \n\
    checkBounds(mouseX, mouseY, d3yR(yLoc) + 15); \n\
  } \n\
 \n\
  var bandXB = [2]; \n\
  bandXB[0] = [MARGINS.left, WIDTH - MARGINS.right]; \n\
  bandXB[1] = [MARGINS.top, HEIGHT - MARGINS.bottom]; \n\
  var heightM = [MARGINS.top + 20, HEIGHT - MARGINS.bottom + 20]; \n\
  //modif Pauline, one more conditon on x \n\
  function checkBounds(x, y, toolYpos) { \n\
    if ( \n\
      bandXB[0][0] <= x && \n\
      bandXB[1][0] <= y && \n\
      x <= bandXB[0][1] && \n\
      y <= bandXB[1][1] && \n\
      heightM[0] <= toolYpos && \n\
      toolYpos <= heightM[1] \n\
    ) { \n\
      $(\"#bandTip\").attr(\"opacity\", 1); \n\
    } else { \n\
      $(\"#bandTip\") \n\
        .attr(\"opacity\", 0) \n\
        .attr(\"x\", \"0px\") \n\
        .attr(\"y\", \"0px\"); \n\
    } \n\
  } \n\
 \n\
  bands \n\
    .append(\"line\") \n\
    .attr(\"class\", \"centerLine\") \n\
    .attr(\"stroke-width\", 1) \n\
    .attr(\"x1\", MARGINS.left) \n\
    .attr(\"y1\", 0) \n\
    .attr(\"x2\", WIDTH - MARGINS.right) \n\
    .attr(\"y2\", 0) \n\
    .attr(\"transform\", function() { \n\
      return \"translate(0,\" + d3yR(0) + \")\"; \n\
    }); \n\
 \n\
  bands \n\
    .append(\"svg:line\") \n\
    .attr(\"class\", \"border\") \n\
    .attr(\"stroke\", \"black\") \n\
    .attr(\"stroke-width\", 1) \n\
    .attr(\"fill\", \"none\") \n\
    .attr(\"x1\", MARGINS.left) \n\
    .attr(\"y1\", MARGINS.top) \n\
    .attr(\"x2\", WIDTH - MARGINS.right) \n\
    .attr(\"y2\", MARGINS.top); \n\
  bands \n\
    .append(\"svg:line\") \n\
    .attr(\"class\", \"border\") \n\
    .attr(\"stroke\", \"black\") \n\
    .attr(\"stroke-width\", 1) \n\
    .attr(\"fill\", \"none\") \n\
    .attr(\"x1\", WIDTH - MARGINS.right) \n\
    .attr(\"y1\", MARGINS.top) \n\
    .attr(\"x2\", WIDTH - MARGINS.right) \n\
    .attr(\"y2\", HEIGHT - MARGINS.top); \n\
  var bandTip = bands.append(\"svg\").attr(\"id\", \"bandTip\"); \n\
  bandTip \n\
    .append(\"circle\") \n\
    .attr(\"cx\", \"6\") \n\
    .attr(\"cy\", \"16\") \n\
    .attr(\"fill\", \"rgba(0, 0, 0, .6)\") \n\
    .attr(\"r\", \"6\"); \n\
  bandTip \n\
    .append(\"polygon\") \n\
    .attr(\"points\", \"15, 16, 22, 10, 22, 22 \") \n\
    .attr(\"fill\", \"#333\"); \n\
  bandTip \n\
    .append(\"rect\") \n\
    .attr(\"width\", 88) \n\
    .attr(\"height\", 28) \n\
    .attr(\"x\", 21) \n\
    .attr(\"y\", 2) \n\
    .attr(\"rx\", 9) \n\
    .attr(\"ry\", 9) \n\
    .style(\"fill\", \"rgba(0,0,0,.8)\") \n\
    .attr(\"stroke-width\", 2) \n\
    .attr(\"stroke\", \"steelblue\"); \n\
  bandTip \n\
    .append(\"text\") \n\
    .attr(\"x\", 65) \n\
    .attr(\"y\", 21) \n\
    .attr(\"fill\", \"white\") \n\
    .attr(\"font-size\", \"14px\") \n\
    .attr(\"text-anchor\", \"middle\") \n\
    .attr(\"id\", \"bandText\"); \n\
  $(\"#bandTip\").attr(\"opacity\", 0); \n\
 \n\
  // Dos plot code \n\
  var dos = d3.select(\"#dos_wrapper\"), \n\
    WIDTH = 300, \n\
    HEIGHT = 550, \n\
    MARGINS = { \n\
      top: 30, \n\
      right: 20, \n\
      bottom: 30, \n\
      left: 10 \n\
    }, \n\
    d3xR2 = d3.scale \n\
      .linear() \n\
      .range([MARGINS.left, WIDTH - MARGINS.right]) \n\
      .domain([dosX[0], dosX[1]]), \n\
    xAxis_dos = d3.svg \n\
      .axis() \n\
      .scale(d3xR2) \n\
      .tickSubdivide(true) \n\
      .ticks(5)  \n\
      .tickSize(5), \n\
    yAxis_dos = d3.svg \n\
      .axis() \n\
      .scale(d3yR) \n\
      .orient(\"left\") \n\
      .tickPadding(10) \n\
      .tickSubdivide(true); \n\
 \n\
  dos.on(\"mousemove\", mousemove).on(\"mouseover\", function() { \n\
    $(\"#bandTip\") \n\
      .attr(\"opacity\", 0) \n\
      .attr(\"x\", 0 + \"px\") \n\
      .attr(\"y\", 0 + \"px\"); \n\
  }); \n\
 \n\
  var bisect = d3.bisector(function(d) { \n\
    return d; \n\
  }).right; \n\
  var dosXB = [2]; \n\
 \n\
  dosXB[0] = [MARGINS.left, WIDTH - MARGINS.right]; \n\
  dosXB[1] = [MARGINS.top, HEIGHT - MARGINS.bottom]; \n\
  var widthM = [MARGINS.left + 20, WIDTH - MARGINS.right + 20]; \n\
 \n\
  function checkBounds2(x, y, toolXpos) { \n\
    if ( \n\
      dosXB[0][0] <= x && \n\
      x <= dosXB[0][1] && \n\
      dosXB[1][0] <= y && \n\
      y <= dosXB[1][1] && \n\
      widthM[0] <= toolXpos && \n\
      toolXpos <= widthM[1] \n\
    ) { \n\
      $(\"#dosTip\").attr(\"opacity\", 1); \n\
      $(\"#dosText\").css(\"opacity\", 1); \n\
    } else { \n\
      $(\"#dosTip\") \n\
        .attr(\"opacity\", 0) \n\
        .attr(\"cx\", 0 + \"px\") \n\
        .attr(\"cy\", 0 + \"px\"); \n\
      $(\"#dosText\").css(\"opacity\", 0); \n\
    } \n\
  } \n\
  dos.on(\"mousemove\", function() { \n\
    var mouseX = d3.mouse(this)[0], \n\
      mouseY = d3.mouse(this)[1]; \n\
 \n\
    var xPos = d3xR2.invert(mouseX), \n\
      yPos = d3yR.invert(mouseY), \n\
      yIndex = Math.min(dosYArray.length - 1, bisect(dosYArray, yPos)), \n\
      dosUnsorted = dosData[yIndex], \n\
      dosSorted = dosUnsorted.slice(0).sort(function(a, b) { \n\
        return a - b; \n\
      }), \n\
      i = Math.max(0, bisect(dosSorted, xPos) - 1), \n\
      d0 = dosSorted[i], \n\
      d1 = dosSorted[Math.min(i + 1, numDos - 1)], \n\
      xIndex = xPos - d0 > d1 - xPos ? Math.min(i + 1, numDos - 1) : i, \n\
      realIndex = dosUnsorted.indexOf(dosSorted[xIndex]), \n\
      name = dosNames[realIndex], \n\
      x1 = dosData[Math.max(yIndex - 1, 0)][realIndex], \n\
      x2 = dosData[yIndex][realIndex], \n\
      y1 = dosYArray[Math.max(yIndex - 1, 0)], \n\
      y2 = dosYArray[yIndex], \n\
      yLoc = y1 + (xPos - x1) * (y2 - y1) / (x2 - x1); // Linear interpolation \n\
 \n\
    $(\"#dosTip\") \n\
      .attr(\"cx\", d3xR2(x2) + \"px\") \n\
      .attr(\"cy\", d3yR(yPos) + \"px\"); \n\
    $(\"#dosText\") \n\
      .html(name + \": \" + x2.toFixed(4) + \" States/eV\") \n\
      .css(\"border-color\", dosColor[realIndex]); \n\
    checkBounds2(mouseX, mouseY, d3xR2(x2) + 21); \n\
  }); \n\
 \n\
  dos \n\
    .append(\"clipPath\") \n\
    .attr(\"id\", \"clip2\") \n\
    .append(\"rect\") \n\
    .attr(\"width\", WIDTH - MARGINS.left - MARGINS.right) \n\
    .attr(\"height\", HEIGHT - MARGINS.top - MARGINS.bottom) \n\
    .attr(\"x\", MARGINS.left) \n\
    .attr(\"y\", MARGINS.top); \n\
 \n\
  dos \n\
    .append(\"svg:g\") \n\
    .attr(\"class\", \"xAxisDos\") \n\
    .attr(\"transform\", \"translate(0,\" + (HEIGHT - MARGINS.bottom) + \")\") \n\
    .call(xAxis_dos); \n\
 \n\
  dos \n\
    .append(\"text\") \n\
    .attr(\"class\", \"titles\") \n\
    .attr(\"x\", WIDTH / 2) \n\
    .attr(\"y\", 15) \n\
    .attr(\"text-anchor\", \"middle\") \n\
    .text(\"Density of States (States/eV)\"); \n\
 \n\
  if (spinP == true) { \n\
    // Add Legend Entries \n\
    for (var i = 0; i < legendNames.length; ++i) { \n\
      $(\"#dosLegend\").append( \n\
        '<text  class=\"legendText\">' + \n\
          legendNames[i] + \n\
          '<line class=\"legendLine\" style=\"border-color:' + \n\
          dosColor[i] + \n\
          ';\"/></text>' \n\
      ); \n\
    } \n\
  } else { \n\
    // Add Legend Entries \n\
    for (var i = 0; i < dosNames.length; ++i) { \n\
      $(\"#dosLegend\").append( \n\
        '<text  class=\"legendText\">' + \n\
          dosNames[i] + \n\
          '<line class=\"legendLine\" style=\"border-color:' + \n\
          dosColor[i] + \n\
          ';\"/></text>' \n\
      ); \n\
    } \n\
  } \n\
  $(\"#dosLegend\").fadeIn(350); \n\
 \n\
  var line2 = d3.svg \n\
    .line() \n\
    .interpolate(\"linear\") \n\
    .x(function(d) { \n\
      return d3xR2(d.x); \n\
    }) \n\
    .y(function(d) { \n\
      return d3yR(d.y); \n\
    }); \n\
  dos \n\
    .selectAll(\".dos_data\") \n\
    .data(dosArray) \n\
    .enter() \n\
    .append(\"path\") \n\
    .attr(\"class\", \"dos_data\") \n\
    .attr(\"stroke\", function(d, i) { \n\
      return dosColor[i]; \n\
    }) \n\
    .attr(\"stroke-width\", 1) \n\
    .attr(\"fill\", \"none\") \n\
    .attr(\"clip-path\", \"url(#clip2)\") \n\
    .attr(\"d\", line2); \n\
 \n\
  dos \n\
    .append(\"line\") \n\
    .attr(\"class\", \"centerLine\") \n\
    .attr(\"stroke-width\", 1) \n\
    .attr(\"x1\", MARGINS.left) \n\
    .attr(\"y1\", 0) \n\
    .attr(\"x2\", WIDTH - MARGINS.right) \n\
    .attr(\"y2\", 0) \n\
    .attr(\"transform\", function() { \n\
      return \"translate(0,\" + d3yR(0) + \")\"; \n\
    }); \n\
 \n\
  dos \n\
    .append(\"svg:line\") \n\
    .attr(\"class\", \"border\") \n\
    .attr(\"stroke\", \"black\") \n\
    .attr(\"stroke-width\", 1) \n\
    .attr(\"fill\", \"none\") \n\
    .attr(\"x1\", MARGINS.left) \n\
    .attr(\"y1\", MARGINS.top) \n\
    .attr(\"x2\", WIDTH - MARGINS.right) \n\
    .attr(\"y2\", MARGINS.top); \n\
  dos \n\
    .append(\"svg:line\") \n\
    .attr(\"class\", \"border\") \n\
    .attr(\"stroke\", \"black\") \n\
    .attr(\"stroke-width\", 1) \n\
    .attr(\"fill\", \"none\") \n\
    .attr(\"x1\", WIDTH - MARGINS.right) \n\
    .attr(\"y1\", MARGINS.top) \n\
    .attr(\"x2\", WIDTH - MARGINS.right) \n\
    .attr(\"y2\", HEIGHT - MARGINS.top); \n\
  dos \n\
    .append(\"svg:line\") \n\
    .attr(\"class\", \"border\") \n\
    .attr(\"stroke\", \"black\") \n\
    .attr(\"stroke-width\", 1) \n\
    .attr(\"fill\", \"none\") \n\
    .attr(\"x1\", MARGINS.left) \n\
    .attr(\"y1\", MARGINS.top) \n\
    .attr(\"x2\", MARGINS.left) \n\
    .attr(\"y2\", HEIGHT - MARGINS.top); \n\
  dos \n\
    .append(\"circle\") \n\
    .attr(\"id\", \"dosTip\") \n\
    .attr(\"fill\", \"rgba(0, 0, 0, .6)\") \n\
    .attr(\"r\", \"6\"); \n\
  $(\"#dosTip\").attr(\"opacity\", \"0\"); \n\
 \n\
  function zoomFunc() { \n\
    //console.log(d3xR(0), d3xR2(0)); \n\
    bands.selectAll(\".band_lines\").attr(\"d\", line); \n\
    bands.select(\".x.axis\").call(xAxis); \n\
    bands.select(\".y.axis\").call(yAxis); \n\
    bands.select(\".centerLine\").attr(\"transform\", function() { \n\
      return \"translate(0,\" + d3yR(0) + \")\"; \n\
    }); \n\
    bands.selectAll(\".lineLabels\").attr(\"transform\", function(d, i) { \n\
      return \"translate(\" + d3xR(lineLocs[i]) + \",0)\"; \n\
    }); \n\
    dos.selectAll(\"path.dos_data\").attr(\"d\", line2); \n\
    dos.select(\".xAxisDos\").call(xAxis_dos); \n\
    dos.select(\".centerLine\").attr(\"transform\", function() { \n\
      return \"translate(0,\" + d3yR(0) + \")\"; \n\
    }); \n\
  } \n\
  var width = WIDTH - MARGINS.left - MARGINS.right, \n\
    height = HEIGHT - MARGINS.top - MARGINS.bottom, \n\
    zoom1b = d3.behavior \n\
      .zoom() \n\
      .x(d3xR) \n\
      .y(d3yR) \n\
      .scaleExtent([1, 100]) \n\
      .on(\"zoom\", zoomFunc), \n\
    zoom1x = d3.behavior \n\
      .zoom() \n\
      .x(d3xR) \n\
      .scaleExtent([1, 100]) \n\
      .on(\"zoom\", zoomFunc), \n\
    zoom1y = d3.behavior \n\
      .zoom() \n\
      .y(d3yR) \n\
      .scaleExtent([1, 100]) \n\
      .on(\"zoom\", zoomFunc), \n\
    zoom2b = d3.behavior \n\
      .zoom() \n\
      .x(d3xR2) \n\
      .y(d3yR) \n\
      .scaleExtent([1, 100]) \n\
      .on(\"zoom\", zoomFunc), \n\
    zoom2x = d3.behavior \n\
      .zoom() \n\
      .x(d3xR2) \n\
      .scaleExtent([1, 100]) \n\
      .on(\"zoom\", zoomFunc), \n\
    zoom2y = d3.behavior \n\
      .zoom() \n\
      .y(d3yR) \n\
      .scaleExtent([1, 100]) \n\
      .on(\"zoom\", zoomFunc), \n\
    zoom1 = [zoom1b, zoom1x, zoom1y], \n\
    zoom2 = [zoom2b, zoom2x, zoom2y]; \n\
 \n\
  // Clicking the Reset Zoom button \n\
  $(\".reset\").click(function() { \n\
    bands.call(zoom1[0].translate([0, 0]).scale(1).event); \n\
    dos.call(zoom2[0].translate([0, 0]).scale(1).event); \n\
  }); \n\
 \n\
  bands.call(zoom1[0]); \n\
  dos.call(zoom2[0]); \n\
 \n\
  $(\"#zoomOptions select\").change(function() { \n\
    var selectV = $(this).val(); \n\
    //console.log(\"test1\"); \n\
    switch (selectV) { \n\
      case \"both\": \n\
        //console.log(\"test2\"); \n\
        bands.call(zoom1[0]); \n\
        dos.call(zoom2[0]); \n\
        break; \n\
      case \"xOnly\": \n\
        bands.call(zoom1[1]); \n\
        dos.call(zoom2[1]); \n\
        break; \n\
      case \"yOnly\": \n\
        bands.call(zoom1[2]); \n\
        dos.call(zoom2[2]); \n\
        break; \n\
    } \n\
  }); \n\
  if (spinP == true) { \n\
    $(\"#spinBandsOptions select\").change(function() { \n\
      var selectVal = $(this).val(); \n\
      switch (selectVal) { \n\
        case \"majority\": \n\
          break; \n\
        case \"minority\": \n\
          break; \n\
      } \n\
    }); \n\
  } \n\
  $(\"#show\").fadeIn(1000); \n\
}); \n\
";
