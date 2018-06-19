// aflowlib_webapp_entry.cpp automatic generated
std::string AFLOW_WEBAPP_ENTRY_JS="\
// aflowlib_webapp_entry.js  \n\
// \n\
// author: Bob Hanson  \n\
// called from aflow_web_interface2.cpp via \n\
// http://aflowlib.duke.edu/users/jmolers/test/entry.php?id=aflow:137d3cb39fa592d3 \n\
 \n\
// BH 5/21/2016 12:57:40 PM \n\
// BH 5/22/2016 10:03:22 PM finalized \n\
// CO 6/23/2017 12:02:22 PM minor changes \n\
 \n\
//var jsmolDir = \"http://aflowlib.duke.edu/users/jmolers/test/jsmol\"; \n\
 \n\
document.write(\"<div id='jmol'></div>\"); \n\
 \n\
//Jmol.db._DirectDatabaseCalls[\"aflowlib.duke.edu\"] = \"%URL\"; // not indicated in JSmolCore.js \n\
//Jmol.db._DirectDatabaseCalls[\"aflow.duke.edu\"] = \"%URL\"; // not indicated in JSmolCore.js \n\
//Jmol.db._DirectDatabaseCalls[\"aflowlib.duke.edu\"] = \"%URL\"; // not indicated in JSmolCore.js \n\
 \n\
if (!self.AFLOW) { \n\
  AFLOW={}; \n\
  AFLOW.version = \"31028\"; \n\
  AFLOW.url_WEB = \"http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/BCC/Ag3Au1Se2_ICSD_171959\"; \n\
  AFLOW.label = \"Ag3Au1Se2_ICSD_171959\"; \n\
  AFLOW.spaceGroupNo = 214; \n\
  AFLOW.spaceGroupName = \"I4_{1}32\"; \n\
  AFLOW.cif_sconv = [10.130,10.130,10.130,90.000,90.000,90.000]; \n\
  AFLOW.cif = [8.773,8.773,8.773,109.471,109.471,109.471]; \n\
  AFLOW.cif_sprim = [8.773,8.773,8.773,109.471,109.471,109.471]; \n\
  AFLOW.baderUnitcell = [-5.065,5.065,5.065,5.065,-5.065,5.065,5.065,5.065,-5.065]; \n\
  AFLOW.baderVSpecies = [\"Ag\",\"Au\",\"Se\"]; \n\
  AFLOW.jsmolDir = \"../test/jsmol\"; \n\
} \n\
AFLOW.url = AFLOW.url_WEB + \"/\" + AFLOW.label; \n\
 \n\
// top left message about the space group \n\
AFLOW.cif2oss = function(params) { \n\
  return \" font echo 14;color echo white; set echo 3% 97%; echo AFLOWLIB consortium (AFLOW v\" + AFLOW.version + \") | entry=\" + AFLOW.label+\"  |  | \" \n\
  + (AFLOW.spaceGroupNo ? \"  Spacegroup = \" + AFLOW.spaceGroupName + \" (#\" + AFLOW.spaceGroupNo + \")   |\" \n\
  : \"\") +  \" a=\" + params[0] + \"\\u212B, b=\" + params[1] + \"\\u212B, c=\" + params[2] + \"\\u212B    |\"  \n\
  + \" \\u03B1=\" + params[3] + \"\\u00B0, \\u03B2=\" + params[4] + \"\\u00B0, \\u03B3=\" + params[5] + \"\\u00B0 ; \"; \n\
} \n\
 \n\
// isosurface creation \n\
AFLOW.colors = [\"red\", \"green\", \"yellow\", \"blue\", \"orange\", \"white\", \"purple\", \"brown\", \"pink\"]; \n\
AFLOW.iso2oss = function(element, cutoff, index) { \n\
  return \"  ISOSURFACE  \" + element + \" '\" + AFLOW.url +\"_Bader_\" + cutoff + \"_\" + element + \".jvxl';\" \n\
    + \"isosurface mesh; color isosurface \" + AFLOW.colors[index] + \" translucent;\" \n\
    + \"VAR charges = load('\" + AFLOW.url + \"_abader.out\" + \"').split('=')[2].split('(')[1].split(','); {*}.label = charges; label %[label];\"; \n\
} \n\
 \n\
// load command, with optional parameters \n\
AFLOW.load = function(root, loadparams) { \n\
  var key = root+loadparams ; \n\
  return \"if (key != '\"+ key + \"') {load '\" + AFLOW.url + root + \".cif' \"  \n\
    + (loadparams ? loadparams : \"packed\") + \";\" \n\
    + AFLOW.cif2oss(AFLOW[\"cif\"+root]) + \";key='\"+key+\"'}\";  \n\
} \n\
 \n\
// JSmol info \n\
AFLOW.Info = { \n\
  width:  Math.floor(window.innerWidth*0.60), \n\
  height: Math.floor(window.innerWidth*0.60), \n\
  debug: false, \n\
  color: \"black\", \n\
  addSelectionOptions: false, \n\
  serverURL: \"https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php\", \n\
  use: \"HTML5\", \n\
  j2sPath: AFLOW.jsmolDir + \"/j2s\", \n\
  readyFunction: function(applet) {Jmol._getElement(applet, \"appletdiv\").style.border=\"1px solid blue\";}, \n\
  script: \"set zoomlarge false; set platformspeed 3;set antialiasDisplay; frank off; set showUnitCellinfo false; \" + AFLOW.load(\"_sconv\"), \n\
  //disableJ2SLoadMonitor: true, \n\
  disableInitialConsole: true \n\
} \n\
 \n\
// thumb wheel widget \n\
AFLOW.input = function(id, val) { \n\
  return \"<input class='dim' id='\" + id + \"' type='number' value='\"+(val||2) + \"' style='width:60px' >\";  \n\
} \n\
 \n\
$(document).ready(function() { \n\
 \n\
  // create header, jmol html, options, and save options, then put them all together at the end \n\
  var header = \"<br /><br /><b>Space Group</b>: \" \n\
     + (AFLOW.spaceGroupNo ? AFLOW.spaceGroupName + \" (#\" + AFLOW.spaceGroupNo + \")\" : \"N/A\") \n\
      \n\
  var jmol = Jmol.getAppletHtml(\"jmolApplet0\", AFLOW.Info) \n\
     \n\
  var options = \"\" \n\
   \n\
  Jmol.setButtonCss(null,\"style='width:180px'\"); \n\
   \n\
  options += \"<b>Relaxed Structure:</b><br /><br />\"      \n\
     + Jmol.jmolButton(jmolApplet0, AFLOW.load(\"\"), \"As calculated\") \n\
           + \"<br />\" \n\
     + Jmol.jmolButton(jmolApplet0, AFLOW.load(\"_sconv\"), \"Standard conventional\") \n\
     + \"<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>\" \n\
           + \"<br />\" \n\
	   + Jmol.jmolButton(jmolApplet0, AFLOW.load(\"_sprim\"), \"Standard primitive\") \n\
     + \"<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>\"; \n\
      \n\
  Jmol.setButtonCss(null,\"style='width:110px'\"); \n\
   \n\
  options += \"<br /><br /><b>Supercell:</b><br /><br />\" \n\
     + Jmol.jmolButton(jmolApplet0, \"key ='';load '' {2 2 2} packed\", \"2x2x2\") \n\
     + Jmol.jmolButton(jmolApplet0, \"key ='';load '' fill 20\", \"20&#8491; box\") \n\
     + \"<p>\" + AFLOW.input(\"dim_1\") + \" X \" + AFLOW.input(\"dim_2\") + \" X \" + AFLOW.input(\"dim_3\") + \"</p>\" \n\
     + \"<input type='button' id='build_button' value='Build' style='width:160px'>\" \n\
     + Jmol.jmolButton(jmolApplet0, AFLOW.load(\"_sconv\") + \";reset;\", \"RESET\"); \n\
 \n\
  Jmol.setButtonCss(null,\"style='width:110px'\")   \n\
 \n\
  options += \"<br /><br /><b>Visualization:</b><br /><br />\" \n\
   \n\
     + Jmol.jmolButton(jmolApplet0, \"spacefill only;spacefill 23%;wireframe 0.15\",\"Ball & Stick\") \n\
     + Jmol.jmolButton(jmolApplet0, \"spacefill #alt:SETTING van der Waals Spheres\", \"Spacefill\") \n\
           + \"<br />\" \n\
     + Jmol.jmolCheckbox(jmolApplet0, \"spin on\",\"spin off\",\"Rotation\") \n\
     + Jmol.jmolCheckbox(jmolApplet0, \"label  %a \",\"labels off\",\"Labels\") \n\
     + Jmol.jmolCheckbox(jmolApplet0, \"background white\",\"background black\", \"Background\"); \n\
     Jmol.setButtonCss(null,\"style='width:30px'\");     \n\
            \n\
     options += \"<br />axis: \"  \n\
     + Jmol.jmolButton(jmolApplet0, \"moveto axis a\",\"a\") \n\
     + Jmol.jmolButton(jmolApplet0, \"moveto axis b\",\"b\") \n\
     + Jmol.jmolButton(jmolApplet0, \"moveto axis c\",\"c\"); \n\
     Jmol.setButtonCss(null,\"style='width:110px'\");     \n\
     options += Jmol.jmolButton(jmolApplet0, \"reset\",\"RESET\"); \n\
 \n\
  options += \"<br /><br /><b>Crystallographic Planes:</b><br />\" \n\
     + \"<p> h:\" + AFLOW.input(\"plane_1\") + \" k:\" + AFLOW.input(\"plane_2\") + \" l:\" + AFLOW.input(\"plane_3\") + \"</p>\" \n\
     + \"<input type='button' id='plane_button' value='Show plane' style='width:160px'>\" \n\
     + Jmol.jmolButton(jmolApplet0, \"isosurface off\", \"RESET\"); \n\
 \n\
  //BEGIN BADER ISOSURFACES \n\
  if (AFLOW.baderUnitcell) { \n\
 \n\
	  options += \"<br /><br /><b>Bader Isosurfaces:</b><br />\"; \n\
 \n\
	  Jmol.setButtonCss(null,\"style='width:50px'\"); \n\
    var cutoffs = [20, 30, 40, 50];  \n\
    var nSpecies = AFLOW.baderVSpecies.length; \n\
    for (i = 0; i < cutoffs.length; i++) { \n\
      var cutoff = cutoffs[i]; \n\
 	    options += \"<br />Cutoff = 0.\" + cutoff + \"<br />\"; \n\
      var unitcell = \"packed UNITCELL[\"+ AFLOW.baderUnitcell.join(\",\") + \"]\";      \n\
  	  for(var j = 0; j < nSpecies; j++) { \n\
        var spec = AFLOW.baderVSpecies[j]; \n\
   	    options += Jmol.jmolButton(jmolApplet0, AFLOW.load(\"\", unitcell) \n\
  	         + \";isosurface delete;\" + AFLOW.iso2oss(spec, cutoff, j), spec); \n\
  	  } \n\
  	  var s = \"isosurface delete;\" + AFLOW.load(\"\", unitcell); \n\
  	  for(var j = 0; j < nSpecies; j++) \n\
  	    s += AFLOW.iso2oss(AFLOW.baderVSpecies[j], cutoff, j); \n\
  	  options += Jmol.jmolButton(jmolApplet0, s, \"All\"); \n\
    } \n\
    options += \"<br />\"; \n\
	} \n\
  //END BADER ISOSURFACES \n\
   \n\
  Jmol.setButtonCss(null,\"style='width:140px'\"); \n\
 \n\
  var saveoptions = \"<b>Save:</b>\" \n\
     + Jmol.jmolButton(jmolApplet0, \"write FILE ?\",\"CIF FILE\") \n\
     + Jmol.jmolButton(jmolApplet0, \"write STATE ?.spt\",\"STATE\") \n\
     + Jmol.jmolButton(jmolApplet0, \"write IMAGE ?.jpg\",\"JPG\") \n\
     + Jmol.jmolButton(jmolApplet0, \"write IMAGE ?.png\",\"PNG\") \n\
     + Jmol.jmolButton(jmolApplet0, \"write PNGJ ?.png\",\"PNG+Jmol\"); \n\
            \n\
  var html =  header + \"<table><tr><td align=center valign=top>\" + jmol + \"</td>\" \n\
            + \"<td><form>\" + options + \"</form></td></tr>\" \n\
            + \"<tr><td align=center>\" + saveoptions + \"</td></tr></table>\"; \n\
             \n\
  $(\"#jmol\").html(html); \n\
   \n\
  // set non-Jmol button click events \n\
  $(\"#build_button\").click(function() { \n\
    var scriptCommand = \"key='';load '' {\" +  $(\"#dim_1\").val() + \" \" + $(\"#dim_2\").val() + \" \" + $(\"#dim_3\").val() + \"} packed;\"; \n\
    Jmol.script(jmolApplet0, scriptCommand); \n\
  }); \n\
  $(\"#plane_button\").click(function() { \n\
      // BH removed: key='';load '' {444 555 -1} packed;  -- why reload?   \n\
    var scriptCommand = \"isosurface ID 'hklplane' hkl {\" +  $(\"#plane_1\").val() + \" \" + $(\"#plane_2\").val() + \" \" + $(\"#plane_3\").val() + \"} colorscheme sets translucent 0.5 green\"; \n\
    Jmol.script(jmolApplet0, scriptCommand); \n\
  }); \n\
}) \n\
";
