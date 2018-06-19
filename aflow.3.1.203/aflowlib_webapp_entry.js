// aflowlib_webapp_entry.js 
//
// author: Bob Hanson 
// called from aflow_web_interface2.cpp via
// http://aflowlib.duke.edu/users/jmolers/test/entry.php?id=aflow:137d3cb39fa592d3

// BH 5/21/2016 12:57:40 PM
// BH 5/22/2016 10:03:22 PM finalized
// CO 6/23/2017 12:02:22 PM minor changes

//var jsmolDir = "http://aflowlib.duke.edu/users/jmolers/test/jsmol";

document.write("<div id='jmol'></div>");

//Jmol.db._DirectDatabaseCalls["aflowlib.duke.edu"] = "%URL"; // not indicated in JSmolCore.js
//Jmol.db._DirectDatabaseCalls["aflow.duke.edu"] = "%URL"; // not indicated in JSmolCore.js
//Jmol.db._DirectDatabaseCalls["aflowlib.duke.edu"] = "%URL"; // not indicated in JSmolCore.js

if (!self.AFLOW) {
  AFLOW={};
  AFLOW.version = "31028";
  AFLOW.url_WEB = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/BCC/Ag3Au1Se2_ICSD_171959";
  AFLOW.label = "Ag3Au1Se2_ICSD_171959";
  AFLOW.spaceGroupNo = 214;
  AFLOW.spaceGroupName = "I4_{1}32";
  AFLOW.cif_sconv = [10.130,10.130,10.130,90.000,90.000,90.000];
  AFLOW.cif = [8.773,8.773,8.773,109.471,109.471,109.471];
  AFLOW.cif_sprim = [8.773,8.773,8.773,109.471,109.471,109.471];
  AFLOW.baderUnitcell = [-5.065,5.065,5.065,5.065,-5.065,5.065,5.065,5.065,-5.065];
  AFLOW.baderVSpecies = ["Ag","Au","Se"];
  AFLOW.jsmolDir = "../test/jsmol";
}
AFLOW.url = AFLOW.url_WEB + "/" + AFLOW.label;

// top left message about the space group
AFLOW.cif2oss = function(params) {
  return " font echo 14;color echo white; set echo 3% 97%; echo AFLOWLIB consortium (AFLOW v" + AFLOW.version + ") | entry=" + AFLOW.label+"  |  | "
  + (AFLOW.spaceGroupNo ? "  Spacegroup = " + AFLOW.spaceGroupName + " (#" + AFLOW.spaceGroupNo + ")   |"
  : "") +  " a=" + params[0] + "\u212B, b=" + params[1] + "\u212B, c=" + params[2] + "\u212B    |" 
  + " \u03B1=" + params[3] + "\u00B0, \u03B2=" + params[4] + "\u00B0, \u03B3=" + params[5] + "\u00B0 ; ";
}

// isosurface creation
AFLOW.colors = ["red", "green", "yellow", "blue", "orange", "white", "purple", "brown", "pink"];
AFLOW.iso2oss = function(element, cutoff, index) {
  return "  ISOSURFACE  " + element + " '" + AFLOW.url +"_Bader_" + cutoff + "_" + element + ".jvxl';"
    + "isosurface mesh; color isosurface " + AFLOW.colors[index] + " translucent;"
    + "VAR charges = load('" + AFLOW.url + "_abader.out" + "').split('=')[2].split('(')[1].split(','); {*}.label = charges; label %[label];";
}

// load command, with optional parameters
AFLOW.load = function(root, loadparams) {
  var key = root+loadparams ;
  return "if (key != '"+ key + "') {load '" + AFLOW.url + root + ".cif' " 
    + (loadparams ? loadparams : "packed") + ";"
    + AFLOW.cif2oss(AFLOW["cif"+root]) + ";key='"+key+"'}"; 
}

// JSmol info
AFLOW.Info = {
  width:  Math.floor(window.innerWidth*0.60),
  height: Math.floor(window.innerWidth*0.60),
  debug: false,
  color: "black",
  addSelectionOptions: false,
  serverURL: "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
  use: "HTML5",
  j2sPath: AFLOW.jsmolDir + "/j2s",
  readyFunction: function(applet) {Jmol._getElement(applet, "appletdiv").style.border="1px solid blue";},
  script: "set zoomlarge false; set platformspeed 3;set antialiasDisplay; frank off; set showUnitCellinfo false; " + AFLOW.load("_sconv"),
  //disableJ2SLoadMonitor: true,
  disableInitialConsole: true
}

// thumb wheel widget
AFLOW.input = function(id, val) {
  return "<input class='dim' id='" + id + "' type='number' value='"+(val||2) + "' style='width:60px' >"; 
}

$(document).ready(function() {

  // create header, jmol html, options, and save options, then put them all together at the end
  var header = "<br /><br /><b>Space Group</b>: "
     + (AFLOW.spaceGroupNo ? AFLOW.spaceGroupName + " (#" + AFLOW.spaceGroupNo + ")" : "N/A")
     
  var jmol = Jmol.getAppletHtml("jmolApplet0", AFLOW.Info)
    
  var options = ""
  
  Jmol.setButtonCss(null,"style='width:180px'");
  
  options += "<b>Relaxed Structure:</b><br /><br />"     
     + Jmol.jmolButton(jmolApplet0, AFLOW.load(""), "As calculated")
           + "<br />"
     + Jmol.jmolButton(jmolApplet0, AFLOW.load("_sconv"), "Standard conventional")
     + "<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>"
           + "<br />"
	   + Jmol.jmolButton(jmolApplet0, AFLOW.load("_sprim"), "Standard primitive")
     + "<a href=http://www.sciencedirect.com/science/article/pii/S0927025610002697>[info]</a>";
     
  Jmol.setButtonCss(null,"style='width:110px'");
  
  options += "<br /><br /><b>Supercell:</b><br /><br />"
     + Jmol.jmolButton(jmolApplet0, "key ='';load '' {2 2 2} packed", "2x2x2")
     + Jmol.jmolButton(jmolApplet0, "key ='';load '' fill 20", "20&#8491; box")
     + "<p>" + AFLOW.input("dim_1") + " X " + AFLOW.input("dim_2") + " X " + AFLOW.input("dim_3") + "</p>"
     + "<input type='button' id='build_button' value='Build' style='width:160px'>"
     + Jmol.jmolButton(jmolApplet0, AFLOW.load("_sconv") + ";reset;", "RESET");

  Jmol.setButtonCss(null,"style='width:110px'")  

  options += "<br /><br /><b>Visualization:</b><br /><br />"
  
     + Jmol.jmolButton(jmolApplet0, "spacefill only;spacefill 23%;wireframe 0.15","Ball & Stick")
     + Jmol.jmolButton(jmolApplet0, "spacefill #alt:SETTING van der Waals Spheres", "Spacefill")
           + "<br />"
     + Jmol.jmolCheckbox(jmolApplet0, "spin on","spin off","Rotation")
     + Jmol.jmolCheckbox(jmolApplet0, "label  %a ","labels off","Labels")
     + Jmol.jmolCheckbox(jmolApplet0, "background white","background black", "Background");
     Jmol.setButtonCss(null,"style='width:30px'");    
           
     options += "<br />axis: " 
     + Jmol.jmolButton(jmolApplet0, "moveto axis a","a")
     + Jmol.jmolButton(jmolApplet0, "moveto axis b","b")
     + Jmol.jmolButton(jmolApplet0, "moveto axis c","c");
     Jmol.setButtonCss(null,"style='width:110px'");    
     options += Jmol.jmolButton(jmolApplet0, "reset","RESET");

  options += "<br /><br /><b>Crystallographic Planes:</b><br />"
     + "<p> h:" + AFLOW.input("plane_1") + " k:" + AFLOW.input("plane_2") + " l:" + AFLOW.input("plane_3") + "</p>"
     + "<input type='button' id='plane_button' value='Show plane' style='width:160px'>"
     + Jmol.jmolButton(jmolApplet0, "isosurface off", "RESET");

  //BEGIN BADER ISOSURFACES
  if (AFLOW.baderUnitcell) {

	  options += "<br /><br /><b>Bader Isosurfaces:</b><br />";

	  Jmol.setButtonCss(null,"style='width:50px'");
    var cutoffs = [20, 30, 40, 50]; 
    var nSpecies = AFLOW.baderVSpecies.length;
    for (i = 0; i < cutoffs.length; i++) {
      var cutoff = cutoffs[i];
 	    options += "<br />Cutoff = 0." + cutoff + "<br />";
      var unitcell = "packed UNITCELL["+ AFLOW.baderUnitcell.join(",") + "]";     
  	  for(var j = 0; j < nSpecies; j++) {
        var spec = AFLOW.baderVSpecies[j];
   	    options += Jmol.jmolButton(jmolApplet0, AFLOW.load("", unitcell)
  	         + ";isosurface delete;" + AFLOW.iso2oss(spec, cutoff, j), spec);
  	  }
  	  var s = "isosurface delete;" + AFLOW.load("", unitcell);
  	  for(var j = 0; j < nSpecies; j++)
  	    s += AFLOW.iso2oss(AFLOW.baderVSpecies[j], cutoff, j);
  	  options += Jmol.jmolButton(jmolApplet0, s, "All");
    }
    options += "<br />";
	}
  //END BADER ISOSURFACES
  
  Jmol.setButtonCss(null,"style='width:140px'");

  var saveoptions = "<b>Save:</b>"
     + Jmol.jmolButton(jmolApplet0, "write FILE ?","CIF FILE")
     + Jmol.jmolButton(jmolApplet0, "write STATE ?.spt","STATE")
     + Jmol.jmolButton(jmolApplet0, "write IMAGE ?.jpg","JPG")
     + Jmol.jmolButton(jmolApplet0, "write IMAGE ?.png","PNG")
     + Jmol.jmolButton(jmolApplet0, "write PNGJ ?.png","PNG+Jmol");
           
  var html =  header + "<table><tr><td align=center valign=top>" + jmol + "</td>"
            + "<td><form>" + options + "</form></td></tr>"
            + "<tr><td align=center>" + saveoptions + "</td></tr></table>";
            
  $("#jmol").html(html);
  
  // set non-Jmol button click events
  $("#build_button").click(function() {
    var scriptCommand = "key='';load '' {" +  $("#dim_1").val() + " " + $("#dim_2").val() + " " + $("#dim_3").val() + "} packed;";
    Jmol.script(jmolApplet0, scriptCommand);
  });
  $("#plane_button").click(function() {
      // BH removed: key='';load '' {444 555 -1} packed;  -- why reload?  
    var scriptCommand = "isosurface ID 'hklplane' hkl {" +  $("#plane_1").val() + " " + $("#plane_2").val() + " " + $("#plane_3").val() + "} colorscheme sets translucent 0.5 green";
    Jmol.script(jmolApplet0, scriptCommand);
  });
})
