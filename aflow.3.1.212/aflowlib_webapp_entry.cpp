// aflowlib_webapp_entry.cpp automatic generated
std::string AFLOW_WEBAPP_ENTRY_JS="\
// aflowlib_webapp_entry.js \n\
// \n\
// author: Bob Hanson \n\
// called from aflow_web_interface2.cpp via \n\
// http://aflowlib.duke.edu/users/jmolers/test/entry.php?id=aflow:137d3cb39fa592d3 \n\
 \n\
// BH 5/21/2016 12:57:40 PM \n\
// BH 5/22/2016 10:03:22 PM finalized \n\
// CO 6/23/2017 12:02:22 PM minor changes \n\
// PC 7/23/2018 Add Symmetry Operations \n\
 \n\
//var jsmolDir = \"http://aflowlib.duke.edu/users/jmolers/test/jsmol\"; \n\
 \n\
document.write(\"<div id='jmol'></div>\"); \n\
 \n\
//Jmol.db._DirectDatabaseCalls[\"aflowlib.duke.edu\"] = \"%URL\"; // not indicated in JSmolCore.js \n\
//Jmol.db._DirectDatabaseCalls[\"aflow.duke.edu\"] = \"%URL\"; // not indicated in JSmolCore.js \n\
//Jmol.db._DirectDatabaseCalls[\"aflowlib.duke.edu\"] = \"%URL\"; // not indicated in JSmolCore.js \n\
 \n\
//function testfunction() { \n\
//  var all = document.getElementById('symop-container').querySelectorAll('*'); \n\
// all.find('input[type=checkbox]:checked').removeAttr('checked'); \n\
//console.log('hello worl!'); \n\
//} \n\
 \n\
if (!self.AFLOW) { \n\
  AFLOW = {}; \n\
  AFLOW.version = \"31028\"; \n\
  AFLOW.url_WEB = \n\
    \"http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/BCC/Ag3Au1Se2_ICSD_171959\"; \n\
  AFLOW.label = \"Ag3Au1Se2_ICSD_171959\"; \n\
  AFLOW.spaceGroupNo = 214; \n\
  AFLOW.spaceGroupName = \"I4_{1}32\"; \n\
  AFLOW.cif_sconv = [10.13, 10.13, 10.13, 90.0, 90.0, 90.0]; \n\
  AFLOW.cif = [8.773, 8.773, 8.773, 109.471, 109.471, 109.471]; \n\
  AFLOW.cif_sprim = [8.773, 8.773, 8.773, 109.471, 109.471, 109.471]; \n\
  AFLOW.baderUnitcell = [ \n\
    -5.065, \n\
    5.065, \n\
    5.065, \n\
    5.065, \n\
    -5.065, \n\
    5.065, \n\
    5.065, \n\
    5.065, \n\
    -5.065 \n\
  ]; \n\
  AFLOW.baderVSpecies = [\"Ag\", \"Au\", \"Se\"]; \n\
  AFLOW.jsmolDir = \"../test/jsmol\"; \n\
} \n\
AFLOW.url = AFLOW.url_WEB + \"/\" + AFLOW.label; \n\
 \n\
// top left message about the space group \n\
AFLOW.cif2oss = function(params) { \n\
  return ( \n\
    \" font echo 14;color echo white; set echo 3% 97%; echo AFLOW.org consortium) (AFLOW v\" + \n\
    AFLOW.version + \n\
    \") | entry=\" + \n\
    AFLOW.label + \n\
    \"  |  | \" + \n\
    (AFLOW.spaceGroupNo \n\
      ? \"  Spacegroup = \" + \n\
        AFLOW.spaceGroupName + \n\
        \" (#\" + \n\
        AFLOW.spaceGroupNo + \n\
        \")   |\" \n\
      : \"\") + \n\
    \" a=\" + \n\
    params[0] + \n\
    \"\\u212B, b=\" + \n\
    params[1] + \n\
    \"\\u212B, c=\" + \n\
    params[2] + \n\
    \"\\u212B    |\" + \n\
    \" \\u03B1=\" + \n\
    params[3] + \n\
    \"\\u00B0, \\u03B2=\" + \n\
    params[4] + \n\
    \"\\u00B0, \\u03B3=\" + \n\
    params[5] + \n\
    \"\\u00B0 ; \" \n\
  ); \n\
}; \n\
 \n\
// isosurface creation \n\
AFLOW.colors = [ \n\
  \"red\", \n\
  \"green\", \n\
  \"yellow\", \n\
  \"blue\", \n\
  \"orange\", \n\
  \"white\", \n\
  \"purple\", \n\
  \"brown\", \n\
  \"pink\" \n\
]; \n\
AFLOW.iso2oss = function(element, cutoff, index) { \n\
  return ( \n\
    \"  ISOSURFACE  \" + \n\
    element + \n\
    \" '\" + \n\
    AFLOW.url + \n\
    \"_Bader_\" + \n\
    cutoff + \n\
    \"_\" + \n\
    element + \n\
    \".jvxl';\" + \n\
    \"isosurface mesh; color isosurface \" + \n\
    AFLOW.colors[index] + \n\
    \" translucent;\" + \n\
    \"VAR charges = load('\" + \n\
    AFLOW.url + \n\
    \"_abader.out\" + \n\
    \"').split('=')[2].split('(')[1].split(','); {*}.label = charges; label %[label];\" \n\
  ); \n\
}; \n\
 \n\
// load command, with optional parameters \n\
// PC -- add a save command to reset the last loaded structure in symmetry part \n\
AFLOW.load = function(root, loadparams) { \n\
  var key = root + loadparams; \n\
  return ( \n\
    \"if (key != '\" + \n\
    key + \n\
    \"') {load '\" + \n\
    AFLOW.url + \n\
    root + \n\
    \".cif' \" + \n\
    (loadparams ? loadparams : \"packed\") + \n\
    \";\" + \n\
    AFLOW.cif2oss(AFLOW[\"cif\" + root]) + \n\
    \";key='\" + \n\
    key + \n\
    \"'}; save state restorestate\" \n\
  ); \n\
}; \n\
 \n\
// JSmol info \n\
AFLOW.Info = { \n\
  width: Math.floor(window.innerWidth * 0.6), \n\
  height: Math.floor(window.innerWidth * 0.6), \n\
  debug: false, \n\
  color: \"black\", \n\
  addSelectionOptions: false, \n\
  serverURL: \"https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php\", \n\
  use: \"HTML5\", \n\
  j2sPath: AFLOW.jsmolDir + \"/j2s\", \n\
  readyFunction: function(applet) { \n\
    Jmol._getElement(applet, \"appletdiv\").style.border = \"1px solid blue\"; \n\
  }, \n\
  script: \n\
    \"set zoomlarge false; set platformspeed 3;set antialiasDisplay; frank off; set showUnitCellinfo false; \" + \n\
    AFLOW.load(\"_sconv\"), \n\
  //disableJ2SLoadMonitor: true, \n\
  disableInitialConsole: true \n\
}; \n\
 \n\
// thumb wheel widget \n\
AFLOW.input = function(id, val) { \n\
  return ( \n\
    \"<input class='dim' id='\" + \n\
    id + \n\
    \"' type='number' value='\" + \n\
    (val || 2) + \n\
    \"' style='width:60px' >\" \n\
  ); \n\
}; \n\
 \n\
//symetry operations -- PC 180723 \n\
// \n\
 \n\
function printDigit(x) { \n\
  return Number.parseFloat(x).toFixed(2); \n\
} \n\
 \n\
function multiplyMinusOne(arr) { \n\
  var newArr = arr.map(function(element) { \n\
    return element * -1; \n\
  }); \n\
  return newArr; \n\
} \n\
 \n\
if (AFLOW.sym2json) { \n\
  var numSymOp = AFLOW.sym2json.length; \n\
  AFLOW.symOp = []; \n\
  var namesSymOp = []; \n\
  var mAxis = []; \n\
  for (var i = 0; i < numSymOp; i++) { \n\
    var axis = AFLOW.sym2json[i][\"axis\"]; \n\
    for (var j = 0; j < axis.length; j++) { \n\
      if (!Number.isInteger(axis[j])) { \n\
        axis[j] = printDigit(axis[j]); \n\
        if (axis[j] == 0.0 || axis[j] == -0.0) { \n\
          axis[j] = 0; \n\
        } \n\
      } \n\
    } \n\
    mAxis[0] = -axis[0]; \n\
    mAxis[1] = -axis[1]; \n\
    mAxis[2] = -axis[2]; \n\
    var angle = AFLOW.sym2json[i][\"angle\"]; \n\
    var checkAngle = \"\"; \n\
    if (angle > 180) { \n\
      checkAngle += \", angle: \" + angle; \n\
    } \n\
    var adname = \"\"; \n\
    var uc = \n\
      \"\" + \n\
      AFLOW.sym2json[i][\"Uc\"][0] + \n\
      \",\" + \n\
      AFLOW.sym2json[i][\"Uc\"][1] + \n\
      \",\" + \n\
      AFLOW.sym2json[i][\"Uc\"][2] + \n\
      \"\"; \n\
    var mUc = \n\
      \"\" + \n\
      multiplyMinusOne(AFLOW.sym2json[i][\"Uc\"][0], -1) + \n\
      \",\" + \n\
      multiplyMinusOne(AFLOW.sym2json[i][\"Uc\"][1], -1) + \n\
      \",\" + \n\
      multiplyMinusOne(AFLOW.sym2json[i][\"Uc\"][2], -1) + \n\
      \"\"; \n\
    var script = \"select 1.1; rotateSelected [[\" + uc + \"]] 80;\"; \n\
    var rotoinversionsScript = \n\
      \"select 1.1; rotateSelected [[\" + \n\
      mUc + \n\
      \"]] 80; invertSelected POINT {0 0 0};\"; \n\
    var drawAxis = \n\
      \"draw axis\" + \n\
      i + \n\
      \" fixed {\" + \n\
      mAxis + \n\
      \"/1} {\" + \n\
      axis + \n\
      \"/1} 200; color $axis\" + \n\
      i + \n\
      \" yellow;\"; \n\
    var removeAxis = \"hide $axis\" + i + \";\"; \n\
    var drawPlane = \n\
      \"isosurface ID plane\" + \n\
      i + \n\
      \" PLANE {\" + \n\
      axis + \n\
      \",0}; color $plane\" + \n\
      i + \n\
      \" yellow translucent 0.4;\"; \n\
    var removePlane = \"hide $plane\" + i + \";\"; \n\
    var reflectionScript = \n\
      \"select 1.1; color atoms translucent 0.8; delay 0.2; select 1.1; color atoms translucent 0;select 1.1; invertSelected PLANE {\" + \n\
      axis + \n\
      \",0};\"; \n\
    var reflectionName = \"reflection,  plane: (\" + axis + \",0)\"; \n\
    var inversionName = \"Inversion -I\"; \n\
    var inversionScript = \"select 1.1; invertSelected POINT {0 0 0};\"; \n\
    var ftau = AFLOW.sym2json[i][\"ftau\"]; \n\
    if (ftau[0] > 0.001 || ftau[1] > 0.001 || ftau[2] > 0.001) { \n\
      script += \"delay 0.5; translateSelected {\" + ftau + \"/1};\"; \n\
      inversionScript += \"delay 0.5; translateSelected {\" + ftau + \"/1};\"; \n\
      reflectionScript += \"delay 0.5; translateSelected {\" + ftau + \"/1};\"; \n\
      rotoinversionsScript += \"delay 0.5; translateSelected {\" + ftau + \"/1};\"; \n\
      adname += \" +translation\"; \n\
      inversionName += \" +translation\"; \n\
      reflectionName += \" +translation\"; \n\
    } \n\
    //numerotate the similar operations \n\
    var name = \n\
      AFLOW.sym2json[i][\"type\"] + \n\
      \" \" + \n\
      AFLOW.sym2json[i][\"Schoenflies\"] + \n\
      adname; \n\
    namesSymOp.push(name); \n\
    var count = 0; \n\
    for (var j = 0; j < namesSymOp.length; j++) { \n\
      if (namesSymOp[j] == name) { \n\
        count += 1; \n\
      } \n\
    } \n\
    name += \" \" + count; \n\
    AFLOW.symOp[i] = { \n\
      name: name, \n\
      reflectionName: reflectionName, \n\
      inversionName: inversionName, \n\
      script: script, \n\
      reflectionScript: reflectionScript, \n\
      inversionScript: inversionScript, \n\
      rotoinversionsScript: rotoinversionsScript, \n\
      drawAxis: drawAxis, \n\
      removeAxis: removeAxis, \n\
      axis: axis, \n\
      mAxis: mAxis, \n\
      angle: angle, \n\
      drawPlane: drawPlane, \n\
      removePlane: removePlane, \n\
      checkAngle: checkAngle \n\
    }; \n\
  } \n\
  var rotations = []; \n\
  var rotoinversions = []; \n\
  for (var i = 0; i < AFLOW.symOp.length; i++) { \n\
    if (AFLOW.symOp[i].name.includes(\"rotation\")) { \n\
      rotations.push(AFLOW.symOp[i].name); \n\
    } else if (AFLOW.symOp[i].name.includes(\"inversion\")) { \n\
      rotoinversions.push(AFLOW.symOp[i].name); \n\
    } \n\
  } \n\
  rotations.sort(); \n\
  rotoinversions.sort(); \n\
  var rotationsSorted = []; \n\
  var rotoinversionsSorted = []; \n\
  var allOperations = []; \n\
  for (var i = 0; i < rotations.length; i++) { \n\
    for (var j = 0; j < AFLOW.symOp.length; j++) { \n\
      if (AFLOW.symOp[j].name == rotations[i]) { \n\
        AFLOW.symOp[j].name = AFLOW.symOp[j].name.substring( \n\
          0, \n\
          AFLOW.symOp[j].name.length - 1 \n\
        ); \n\
        rotationsSorted.push(AFLOW.symOp[j]); \n\
        allOperations.push(AFLOW.symOp[j]); \n\
      } \n\
    } \n\
  } \n\
  for (var i = 0; i < rotoinversions.length; i++) { \n\
    for (var j = 0; j < AFLOW.symOp.length; j++) { \n\
      if (AFLOW.symOp[j].name == rotoinversions[i]) { \n\
        AFLOW.symOp[j].name = AFLOW.symOp[j].name.substring( \n\
          0, \n\
          AFLOW.symOp[j].name.length - 1 \n\
        ); \n\
        rotoinversionsSorted.push(AFLOW.symOp[j]); \n\
        allOperations.push(AFLOW.symOp[j]); \n\
      } \n\
    } \n\
  } \n\
} \n\
 \n\
$(document).ready(function() { \n\
  // create header, jmol html, options, and save options, then put them all together at the end \n\
  var header = \n\
    \"<br /><br /><b>Space Group</b>: \" + \n\
    (AFLOW.spaceGroupNo \n\
      ? AFLOW.spaceGroupName + \" (#\" + AFLOW.spaceGroupNo + \")\" \n\
      : \"N/A\"); \n\
 \n\
  var jmol = Jmol.getAppletHtml(\"jmolApplet0\", AFLOW.Info); \n\
 \n\
  var options = \"\"; \n\
 \n\
  Jmol.setButtonCss(null, \"style='width:180px'\"); \n\
 \n\
  options += \n\
    \"<b>Relaxed Structure:</b><br /><br />\" + \n\
    Jmol.jmolButton(jmolApplet0, AFLOW.load(\"\"), \"As calculated\") + \n\
    \"<br />\" + \n\
    Jmol.jmolButton( \n\
      jmolApplet0, \n\
      AFLOW.load(\"_sconv\"), \n\
      \"Standard conventional\" \n\
    ) + \n\
    \"<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>\" + \n\
    \"<br />\" + \n\
    Jmol.jmolButton(jmolApplet0, AFLOW.load(\"_sprim\"), \"Standard primitive\") + \n\
    \"<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>\"; \n\
 \n\
  Jmol.setButtonCss(null, \"style='width:110px'\"); \n\
 \n\
  options += \n\
    \"<br /><br /><b>Supercell:</b><br /><br />\" + \n\
    Jmol.jmolButton(jmolApplet0, \"key ='';load '' {2 2 2} packed\", \"2x2x2\") + \n\
    Jmol.jmolButton(jmolApplet0, \"key ='';load '' fill 20\", \"20&#8491; box\") + \n\
    \"<p>\" + \n\
    AFLOW.input(\"dim_1\") + \n\
    \" X \" + \n\
    AFLOW.input(\"dim_2\") + \n\
    \" X \" + \n\
    AFLOW.input(\"dim_3\") + \n\
    \"</p>\" + \n\
    \"<input type='button' id='build_button' value='Build' style='width:160px'>\" + \n\
    Jmol.jmolButton(jmolApplet0, AFLOW.load(\"_sconv\") + \";reset;\", \"RESET\"); \n\
 \n\
  Jmol.setButtonCss(null, \"style='width:110px'\"); \n\
 \n\
  options += \n\
    \"<br /><br /><b>Visualization:</b><br /><br />\" + \n\
    Jmol.jmolButton( \n\
      jmolApplet0, \n\
      \"spacefill only;spacefill 23%;wireframe 0.15\", \n\
      \"Ball & Stick\" \n\
    ) + \n\
    Jmol.jmolButton( \n\
      jmolApplet0, \n\
      \"spacefill #alt:SETTING van der Waals Spheres\", \n\
      \"Spacefill\" \n\
    ) + \n\
    \"<br />\" + \n\
    Jmol.jmolCheckbox(jmolApplet0, \"spin on\", \"spin off\", \"Rotation\") + \n\
    Jmol.jmolCheckbox(jmolApplet0, \"label  %a \", \"labels off\", \"Labels\") + \n\
    Jmol.jmolCheckbox( \n\
      jmolApplet0, \n\
      \"background white\", \n\
      \"background black\", \n\
      \"Background\" \n\
    ); \n\
  Jmol.setButtonCss(null, \"style='width:30px'\"); \n\
 \n\
  options += \n\
    \"<br />axis: \" + \n\
    Jmol.jmolButton(jmolApplet0, \"moveto axis a\", \"a\") + \n\
    Jmol.jmolButton(jmolApplet0, \"moveto axis b\", \"b\") + \n\
    Jmol.jmolButton(jmolApplet0, \"moveto axis c\", \"c\"); \n\
  Jmol.setButtonCss(null, \"style='width:110px'\"); \n\
  options += Jmol.jmolButton(jmolApplet0, \"reset\", \"RESET\"); \n\
 \n\
  options += \n\
    \"<br /><br /><b>Crystallographic Planes:</b><br />\" + \n\
    \"<p> h:\" + \n\
    AFLOW.input(\"plane_1\") + \n\
    \" k:\" + \n\
    AFLOW.input(\"plane_2\") + \n\
    \" l:\" + \n\
    AFLOW.input(\"plane_3\") + \n\
    \"</p>\" + \n\
    \"<input type='button' id='plane_button' value='Show plane' style='width:160px'>\" + \n\
    Jmol.jmolButton(jmolApplet0, \"isosurface off\", \"RESET\"); \n\
 \n\
  //BEGIN BADER ISOSURFACES \n\
  var bader = \"\"; \n\
  if (AFLOW.baderUnitcell) { \n\
    bader += \"<br /><br /><b>Bader Isosurfaces:</b><br />\"; \n\
 \n\
    Jmol.setButtonCss(null, \"style='width:50px'\"); \n\
    var cutoffs = [20, 30, 40, 50]; \n\
    var nSpecies = AFLOW.baderVSpecies.length; \n\
    for (i = 0; i < cutoffs.length; i++) { \n\
      var cutoff = cutoffs[i]; \n\
      bader += \"<br />Cutoff = 0.\" + cutoff + \"<br />\"; \n\
      var unitcell = \"packed UNITCELL[\" + AFLOW.baderUnitcell.join(\",\") + \"]\"; \n\
      for (var j = 0; j < nSpecies; j++) { \n\
        var spec = AFLOW.baderVSpecies[j]; \n\
        bader += Jmol.jmolButton( \n\
          jmolApplet0, \n\
          AFLOW.load(\"\", unitcell) + \n\
            \";isosurface delete;\" + \n\
            AFLOW.iso2oss(spec, cutoff, j), \n\
          spec \n\
        ); \n\
      } \n\
      var s = \"isosurface delete;\" + AFLOW.load(\"\", unitcell); \n\
      for (var j = 0; j < nSpecies; j++) \n\
        s += AFLOW.iso2oss(AFLOW.baderVSpecies[j], cutoff, j); \n\
      bader += Jmol.jmolButton(jmolApplet0, s, \"All\"); \n\
    } \n\
    bader += \"<br />\"; \n\
  } \n\
  //END BADER ISOSURFACES \n\
 \n\
  //BEGIN SYMMETRY -- PC 180723 \n\
  // \n\
  var symProperties = \"\"; \n\
  if (AFLOW.sym2json) { \n\
    Jmol.setButtonCss(null, \"style='width:200px'\"); \n\
    symProperties += \"<br /><br /><b>Symmetry Operations: </b><br />\"; \n\
    Jmol.setButtonCss(null, \"style='width:60px'\"); \n\
    symProperties += \n\
      \"<p>\" + \n\
      numSymOp + \n\
      \" operations calculated with respect to the 'As Calculated' material.</p>\"; \n\
    symProperties += Jmol.jmolCheckbox( \n\
      jmolApplet0, \n\
      \"load append '' {444 666 -1}; select 2.1; color translucent 0.9 [224 224 224];  frame all; center; \", \n\
      \"zap 2.1;\", \n\
      \"Background supercell (for better visualization)\" \n\
    ); \n\
    symProperties += \"<div id='symop-container' class='symop-container'>\"; \n\
    symProperties += \"<span style='font-size:small'> Identity\"; \n\
    symProperties += \"</span><br />\"; \n\
    for (var k = 0; k < rotationsSorted.length; k++) { \n\
      symProperties += \n\
        \"<span style='font-size:small'>\" + \n\
        rotationsSorted[k].name + \n\
        rotationsSorted[k].checkAngle + \n\
        \",   axis: (\" + \n\
        rotationsSorted[k].axis + \n\
        \"), \"; \n\
      symProperties += Jmol.jmolCheckbox( \n\
        jmolApplet0, \n\
        rotationsSorted[k].drawAxis, \n\
        rotationsSorted[k].removeAxis, \n\
        \"axis\" \n\
      ); \n\
      symProperties += Jmol.jmolButton( \n\
        jmolApplet0, \n\
        rotationsSorted[k].script, \n\
        \"apply\" \n\
      ); \n\
      symProperties += \"</span><br />\"; \n\
    } \n\
    for (var k = 0; k < rotoinversionsSorted.length; k++) { \n\
      if (rotoinversionsSorted[k].name.includes(\"-I\")) { \n\
        symProperties += \"<span style='font-size:small'>\"; \n\
        symProperties += rotoinversionsSorted[k].inversionName; \n\
        symProperties += Jmol.jmolCheckbox( \n\
          jmolApplet0, \n\
          \"polyhedra ID pol {0 0 0} TO [{0 0 0.25/1}, {-0.25 0 0/1}, {0 0.25 0/1}, {0 0 -0.25/1}, {0.25 0 0/1}, {0 -0.25 0/1}] ; color $pol yellow;\", \n\
          \"delete $pol\", \n\
          \"point\" \n\
        ); \n\
        symProperties += Jmol.jmolButton( \n\
          jmolApplet0, \n\
          rotoinversionsSorted[k].inversionScript, \n\
          \"apply\" \n\
        ); \n\
        symProperties += \"</span><br />\"; \n\
      } else if (rotoinversionsSorted[k].name.includes(\"S2\")) { \n\
        symProperties += \"<span style='font-size:small'>\"; \n\
        symProperties += rotoinversionsSorted[k].reflectionName; \n\
        symProperties += Jmol.jmolCheckbox( \n\
          jmolApplet0, \n\
          rotoinversionsSorted[k].drawPlane, \n\
          rotoinversionsSorted[k].removePlane, \n\
          \"plane\" \n\
        ); \n\
        symProperties += Jmol.jmolButton( \n\
          jmolApplet0, \n\
          rotoinversionsSorted[k].reflectionScript, \n\
          \"apply\" \n\
        ); \n\
        symProperties += \"</span><br />\"; \n\
      } else if (rotoinversionsSorted[k].name.includes(\"rotoinversion s\")) { \n\
        symProperties += \"<span style='font-size:small'>\"; \n\
        symProperties += rotoinversionsSorted[k].reflectionName; \n\
        symProperties += Jmol.jmolCheckbox( \n\
          jmolApplet0, \n\
          rotoinversionsSorted[k].drawPlane, \n\
          rotoinversionsSorted[k].removePlane, \n\
          \"plane\" \n\
        ); \n\
        symProperties += Jmol.jmolButton( \n\
          jmolApplet0, \n\
          rotoinversionsSorted[k].reflectionScript, \n\
          \"apply\" \n\
        ); \n\
        symProperties += \"</span><br />\"; \n\
      } \n\
    } \n\
    for (var k = 0; k < rotoinversionsSorted.length; k++) { \n\
      if ( \n\
        rotoinversionsSorted[k].name.includes(\"S3\") || \n\
        rotoinversionsSorted[k].name.includes(\"S4\") || \n\
        rotoinversionsSorted[k].name.includes(\"S6\") \n\
      ) { \n\
        symProperties += \"<span style='font-size:small'>\"; \n\
        symProperties += rotoinversionsSorted[k].name; \n\
        symProperties += \n\
          rotoinversionsSorted[k].checkAngle + \n\
          \",   axis: (\" + \n\
          rotoinversionsSorted[k].axis + \n\
          \"), \"; \n\
        symProperties += Jmol.jmolCheckbox( \n\
          jmolApplet0, \n\
          rotoinversionsSorted[k].drawAxis, \n\
          rotoinversionsSorted[k].removeAxis, \n\
          \"axis\" \n\
        ); \n\
        symProperties += Jmol.jmolCheckbox( \n\
          jmolApplet0, \n\
          \"polyhedra ID pol {0 0 0} TO [{0 0 0.25/1}, {-0.25 0 0/1}, {0 0.25 0/1}, {0 0 -0.25/1}, {0.25 0 0/1}, {0 -0.25 0/1}] ; color $pol yellow;\", \n\
          \"delete $pol\", \n\
          \"point\" \n\
        ); \n\
        symProperties += Jmol.jmolButton( \n\
          jmolApplet0, \n\
          rotoinversionsSorted[k].rotoinversionsScript, \n\
          \"apply\" \n\
        ); \n\
        symProperties += \"</span><br />\"; \n\
      } \n\
    } \n\
    symProperties += \"</div>\"; \n\
    Jmol.setButtonCss(null, \"style='width:110px'\"); \n\
    symProperties += \"<form id='test'>\"; \n\
    symProperties += Jmol.jmolButton( \n\
      jmolApplet0, \n\
      \"restore state restorestate;\", \n\
      \"RESET\", \n\
      (id = \"symop-reset\") \n\
    ); \n\
  } \n\
  // END SYMMETRY -- PC 180723 \n\
 \n\
  Jmol.setButtonCss(null, \"style='width:140px'\"); \n\
 \n\
  var saveoptions = \n\
    \"<b>Save:</b>\" + \n\
    Jmol.jmolButton(jmolApplet0, \"write FILE ?\", \"CIF FILE\") + \n\
    Jmol.jmolButton(jmolApplet0, \"write STATE ?.spt\", \"STATE\") + \n\
    Jmol.jmolButton(jmolApplet0, \"write IMAGE ?.jpg\", \"JPG\") + \n\
    Jmol.jmolButton(jmolApplet0, \"write IMAGE ?.png\", \"PNG\") + \n\
    Jmol.jmolButton(jmolApplet0, \"write PNGJ ?.png\", \"PNG+Jmol\"); \n\
 \n\
  var html = \n\
    header + \n\
    \"<table><tr><td align=center valign=top>\" + \n\
    jmol + \n\
    \"</td>\" + \n\
    \"<td><form>\" + \n\
    options + \n\
    symProperties + \n\
    bader + \n\
    \"</form></td></tr>\" + \n\
    \"<tr><td align=center>\" + \n\
    saveoptions + \n\
    \"</td></tr></table>\"; \n\
 \n\
  $(\"#jmol\").html(html); \n\
 \n\
  var resetSym = document.getElementById(\"symop-reset\"); \n\
  resetSym.addEventListener(\"click\", function() { \n\
    console.log(\"enter function\"); \n\
    var checkboxes = document.getElementsByTagName(\"input\"); \n\
    console.log(\"checkboxes\"); \n\
    for (var i = 0; i < checkboxes.length; i++) { \n\
      checkboxes[i].checked = false; \n\
    } \n\
  }); \n\
 \n\
  // set non-Jmol button click events \n\
  $(\"#build_button\").click(function() { \n\
    var scriptCommand = \n\
      \"key='';load '' {\" + \n\
      $(\"#dim_1\").val() + \n\
      \" \" + \n\
      $(\"#dim_2\").val() + \n\
      \" \" + \n\
      $(\"#dim_3\").val() + \n\
      \"} packed;\"; \n\
    Jmol.script(jmolApplet0, scriptCommand); \n\
  }); \n\
  $(\"#plane_button\").click(function() { \n\
    // BH removed: key='';load '' {444 555 -1} packed;  -- why reload? \n\
    var scriptCommand = \n\
      \"isosurface ID 'hklplane' hkl {\" + \n\
      $(\"#plane_1\").val() + \n\
      \" \" + \n\
      $(\"#plane_2\").val() + \n\
      \" \" + \n\
      $(\"#plane_3\").val() + \n\
      \"} colorscheme sets translucent 0.5 green\"; \n\
    Jmol.script(jmolApplet0, scriptCommand); \n\
  }); \n\
}); \n\
";
