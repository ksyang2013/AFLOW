// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo - Duke

#ifndef _AFLOWLIB_WEB_OUTREACH_CPP_
#define _AFLOWLIB_WEB_OUTREACH_CPP_
#include "aflow.h"

#define TOC 0

#define WEB_PDF  string("http://"+XHOST.AFLOW_MATERIALS_SERVER+"/auro/AUROARTICULA/")
#define WEB_DOI  string("https://doi.org/")

#define THRUST_RECENT_ARTICLES  20
#define THRUST_RECENT_YEARS     5
#define AUTHOR_RECENT_ARTICLES  200
#define AUTHOR_RECENT_YEARS     30

#define MAX_YEAR_PRESENTATIONS  2019

// ******************************************************************************************************************************************************
// _OUTREACH CLASS
// namespace web

_outreach::_outreach() {
  newflag=FALSE;
  wnumber=0;
  year=0;
  vauthor.clear();
  title="";
  journal="";
  arxiv="";
  supplementary="";
  bibtex="";
  place="";
  date="";
  link="";
  type="";
  _isinvited=FALSE;
  host="";
  abstract="";
  pdf="";
  doi="";
  vextra_html.clear();
  vextra_latex.clear();
  vkeyword.clear();
  vsponsor.clear();
  valloy.clear();
}

// destructor
_outreach::~_outreach() {
  free();
}

void _outreach::free() {
}

void _outreach::copy(const _outreach& b) {
  // const _outreach& _outreach::operator=(const _outreach& b) {       // operator=
  free();
  newflag=b.newflag;
  wnumber=b.wnumber;
  year=b.year;
  vauthor.clear();for(uint i=0;i<b.vauthor.size();i++) vauthor.push_back(b.vauthor.at(i));
  title=b.title;
  journal=b.journal;
  arxiv=b.arxiv;
  supplementary=b.supplementary;
  bibtex=b.bibtex;
  place=b.place;
  date=b.date;
  link=b.link;
  type=b.type;
  _isinvited=b._isinvited;
  host=b.host;
  abstract=b.abstract;
  pdf=b.pdf;
  doi=b.doi;
  vextra_html.clear();for(uint i=0;i<b.vextra_html.size();i++) vextra_html.push_back(b.vextra_html.at(i));
  vextra_latex.clear();for(uint i=0;i<b.vextra_latex.size();i++) vextra_latex.push_back(b.vextra_latex.at(i));
  vkeyword.clear();for(uint i=0;i<b.vkeyword.size();i++) vkeyword.push_back(b.vkeyword.at(i));
  vsponsor.clear();for(uint i=0;i<b.vsponsor.size();i++) vsponsor.push_back(b.vsponsor.at(i));
  valloy.clear();for(uint i=0;i<b.valloy.size();i++) valloy.push_back(b.valloy.at(i));
}

// copy
_outreach::_outreach(const _outreach& b) {
  // [OBSOLETE]   free();
  // [OBSOLETE]  *this=b;
  copy(b);
}

// copies xtructures: b=a
const _outreach& _outreach::operator=(const _outreach& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _outreach::clear() {
  _outreach outreach_temp;
  copy(outreach_temp);
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// ARTICLES/PRESENTATIONS
// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

string fixonelabel(const vector<string>& vlabel,string label) {
  for(uint i=0;i<vlabel.size();i+=2)
    if(label==vlabel.at(i)) return vlabel.at(i+1);
  return label;
}

uint fixlabel(const vector<string>& vlabel,vector<string>& vlabelfixed) {
  for(uint i=0;i<vlabelfixed.size();i++) {
    for(uint j=0;j<vlabel.size();j+=2) 
      if(vlabelfixed.at(i)==vlabel.at(j))
	aurostd::StringSubst(vlabelfixed.at(i),vlabelfixed.at(i),vlabel.at(j+1));
    aurostd::StringSubst(vlabelfixed.at(i),"_WEB_PDF_",WEB_PDF);
    aurostd::StringSubst(vlabelfixed.at(i),"_WEB_DOI_",WEB_DOI);
  }
  return vlabelfixed.size();
}

ostream& operator<<(ostream& oss,const _outreach& outreach) {
  // ******************************************************************************************************************************************************
  // ARTICLE
  if(aurostd::substring2bool(outreach.type,"ARTICLE")) {
    bool compact=TRUE;stringstream newline;newline << endl;

    string authors;
    for(uint iauth=0;iauth<outreach.vauthor.size();iauth++) {
      authors+=outreach.vauthor.at(iauth);
      if(outreach.vauthor.size()==2 && iauth==outreach.vauthor.size()-2) authors+=" and ";
      if(outreach.vauthor.size()!=2 && iauth==outreach.vauthor.size()-2) authors+=", and ";
      if(iauth!=outreach.vauthor.size()-2 && iauth!=outreach.vauthor.size()-1) authors+=", ";
    }
    string wnumber=aurostd::utype2string(outreach.wnumber);
    if(wnumber.length()<2) wnumber="0"+wnumber;

    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      // 3rd line AUTHORS
      aurostd::StringSubst(authors,"Csányi","Csanyi");
      oss << aurostd::html2txt(authors) << ", ";
      // 4th line TITLE
      oss << "\"" << aurostd::html2txt(outreach.title) << "\", ";
      // 5th line JOURNAL with year
      oss << "" << aurostd::html2txt(outreach.journal) << ". ";
      //   oss  << endl;
      // EXTRA STUFF FOR PHP WEB
      // 7th line DOI
      if(outreach.doi.length()>0)
	oss << "doi: " << outreach.doi;
      //     // 6th line PDF
      //     if(outreach.pdf.length()>0)
      //       oss << "[<a href="+WEB_PDF+outreach.pdf << "><b>pdf</b></a>] " << endl;
      //     // 6th line SUPPLEMENTARY
      //     if(outreach.supplementary.length()>0)
      //       oss << "[<a href="+WEB_PDF+outreach.supplementary << "><b>suppl</b></a>] " << endl;
      //     // 6th line ARXIV
      //     if(outreach.arxiv.length()>0)
      //       oss << "[<a href="+outreach.arxiv << "><b>arxiv</b></a>] " << endl;
      //     // 6th line LINK
      //     if(outreach.link.length()>0)
      //       oss << "[<a href="+outreach.link << "><b>link</b></a>] " << endl;
      //     // Xth line EXTRA
      //     for(uint ix=0;ix<outreach.vextra_html.size();ix++)
      //       oss << outreach.vextra_html.at(ix) << endl;
      //     oss << "<br>";
      oss << endl;
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      oss << "<li>";
      if(outreach.newflag) {
	if(XHOST.vflag_control.flag("PRINT_MODE::NEW") && outreach.newflag) oss << "<blink><b>NEW </b></blink>";
	if(XHOST.vflag_control.flag("PRINT_MODE::DATE") && outreach.newflag && outreach.date.length()>0) oss << "<b>" << outreach.date << "</b>";
	oss << " "; //endl;
      }
      if(XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
	if(outreach.wnumber==0) {
	  oss << "<span class=\"pubNumber\">" << "Chapter" << ". </span>";
	} else {
	  oss << "<span class=\"pubNumber\">" << outreach.wnumber << ". </span>";
	}
      }
      // 3rd line AUTHORS
      aurostd::StringSubst(authors,"Csányi","Csanyi");
      oss << "<span class=\"pubAuthors\">" <<  authors << "</span>, " << (compact?" ":newline.str());
      // 4th line TITLE
      oss << "<span class=\"pubTitle\">" << outreach.title << "</span>, " << (compact?" ":newline.str());
      // 5th line JOURNAL with year
      oss << "<span class=\"pubJournal\">" << outreach.journal << "</span>. " << (compact?" ":newline.str());
      // EXTRA STUFF FOR PHP WEB
      // 7th line DOI
      if(XHOST.vflag_control.flag("PRINT_MODE::DOI") && outreach.doi.length()>0)
	oss << "[<a href="+WEB_DOI+outreach.doi << "><b>doi" << "=" << outreach.doi << "</b></a>] " << (compact?" ":newline.str());
      // 6th line PDF    
      if(XHOST.vflag_control.flag("PRINT_MODE::PDF") && outreach.pdf.length()>0)
	oss << "[<a href="+WEB_PDF+outreach.pdf << "><b>pdf</b></a>] " << (compact?" ":newline.str());
      // 6th line SUPPLEMENTARY
      if(XHOST.vflag_control.flag("PRINT_MODE::PDF") && outreach.supplementary.length()>0)
	oss << "[<a href="+WEB_PDF+outreach.supplementary << "><b>suppl</b></a>] " << (compact?" ":newline.str());
      // 6th line ARXIV
      if(outreach.arxiv.length()>0)
	oss << "[<a href="+outreach.arxiv << "><b>arxiv</b></a>] " << (compact?" ":newline.str());
      // 6th line LINK
      if(outreach.link.length()>0)
	oss << "[<a href="+outreach.link << "><b>link</b></a>] " << (compact?" ":newline.str());
      // Xth line EXTRA
      for(uint ix=0;ix<outreach.vextra_html.size();ix++) {
	if(!aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_PDF) && !aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_DOI)) {
	  oss << outreach.vextra_html.at(ix) << (compact?" ":newline.str());
	} else {
	  if(aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_DOI) && XHOST.vflag_control.flag("PRINT_MODE::DOI")) {
	    oss << outreach.vextra_html.at(ix) << (compact?" ":newline.str());
	  } else {
	    if(aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_PDF) && XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
	      oss << outreach.vextra_html.at(ix) << (compact?" ":newline.str());
	    }
	  }
	}
      }
      //    oss << "<br>";
      oss << "</li>";
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
      // 3rd line AUTHORS
      authors=aurostd::html2latex(authors)+", ";
      oss << " " << authors;
      // 4th line TITLE
      string title="{\\it "+aurostd::html2latex(outreach.title)+"}, ";
      oss << " " << title;
      // 5th line JOURNAL with year
      string journal=aurostd::html2latex(outreach.journal)+". ";
      oss << " " << journal;
      // EXTRA STUFF FOR LATEX
      bool link=FALSE;
      // 7th line DOI
      if(XHOST.vflag_control.flag("PRINT_MODE::DOI") && outreach.doi.length()>0) {
	string doi="";
	doi="\\ifthenelse{\\equal{\\hyperlinks}{true}}{";
	doi+="{\\newline \\sf \\href{https://doi.org/"+outreach.doi+"}{DOI: "+aurostd::html2latex(outreach.doi)+"}}";
	doi+="}{{\\newline \\sf DOI: "+aurostd::html2latex(outreach.doi)+"}}";
	oss << " " << doi;
	link=TRUE;
      } else {
	if(outreach.wnumber!=0 && outreach.wnumber!=27 && outreach.wnumber!=14 &&
	   outreach.wnumber!=11 && outreach.wnumber!=8 && outreach.wnumber!=2 &&
	   outreach.wnumber!=1) {
	  string doi=""; // ="{\\newline \\sf \\href{https://doi.org/}{DOI: N/A}}";
	  oss << " " << doi;
	  link=TRUE;
	}
      }
      // 7th line PDF
      if((XHOST.vflag_control.flag("PRINT_MODE::PDF") || (!XHOST.vflag_control.flag("PRINT_MODE::PDF") && XHOST.vflag_control.flag("PRINT_MODE::DOI"))) && outreach.doi.length()==0 && outreach.pdf.length()>0) {
	string pdf="";
	pdf="\\ifthenelse{\\equal{\\hyperlinks}{true}}{{\\sf [\\href{"+WEB_PDF+outreach.pdf+"}{pdf}]}}{}";
	oss << " " << pdf;
	link=TRUE;
      }
      if(link==FALSE) {;};//cerr  << wnumber << endl;
      // LATEX EXTRA
      for(uint ix=0;ix<outreach.vextra_latex.size();ix++) {
	if(!aurostd::substring2bool(outreach.vextra_latex.at(ix),WEB_PDF) && !aurostd::substring2bool(outreach.vextra_latex.at(ix),WEB_DOI)) {
	  oss << " " << outreach.vextra_latex.at(ix);
	} else {
	  if(aurostd::substring2bool(outreach.vextra_latex.at(ix),WEB_DOI) && XHOST.vflag_control.flag("PRINT_MODE::DOI")) {
	    oss << " " << outreach.vextra_latex.at(ix);
	  } else {
	    if(aurostd::substring2bool(outreach.vextra_latex.at(ix),WEB_PDF) && XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
	      oss << " " <<  outreach.vextra_latex.at(ix);
	    }
	  }
	}
      }
      oss << "% ART" << wnumber << " - aflow " << string(AFLOW_VERSION) << endl;
    }
  }
  // ***************************************************************************
  // PRESENTATION
  if(aurostd::substring2bool(outreach.type,"PRESENTATION")) {
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      if(outreach._isinvited && outreach.type=="PRESENTATION_TALK") oss << "Invited talk:";
      if(outreach._isinvited && outreach.type=="PRESENTATION_SEMINAR") oss << "Invited seminar:";
      if(outreach._isinvited && outreach.type=="PRESENTATION_COLLOQUIUM") oss << "Invited colloquium:";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PLENARY") oss << "Plenary Speaker:";
      if(outreach._isinvited && outreach.type=="PRESENTATION_KEYNOTE") oss << "Keynote Speaker:";
      if(outreach._isinvited && outreach.type=="PRESENTATION_TUTORIAL") oss << "Tutorial:";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PANELIST") oss << "Invited panelist:";
      oss << " " << aurostd::latex2txt(outreach.title) << "; ";
      oss << "" << aurostd::latex2txt(outreach.place) << ", ";
      oss << "" << aurostd::latex2txt(outreach.date) << ". ";
      //    oss << "" << endl;
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      if(outreach._isinvited && outreach.type=="PRESENTATION_TALK") oss << "<b> Invited talk:</b>";
      if(outreach._isinvited && outreach.type=="PRESENTATION_SEMINAR") oss << "<b> Invited seminar:</b>";
      if(outreach._isinvited && outreach.type=="PRESENTATION_COLLOQUIUM") oss << "<b> Invited colloquium:</b>";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PLENARY") oss << "<b> Plenary Speaker:</b>";
      if(outreach._isinvited && outreach.type=="PRESENTATION_KEYNOTE") oss << "<b> Keynote Speaker:</b>";
      if(outreach._isinvited && outreach.type=="PRESENTATION_TUTORIAL") oss << "<b> Tutorial:</b>";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PANELIST") oss << "<b> Invited panelist:</b>";
      //    oss << "<br>" << endl;
      oss << "<i> " << aurostd::latex2html(outreach.title) << "</i>; " << endl;
      oss << "" << aurostd::latex2html(outreach.place) << ", " << endl;
      oss << "" << aurostd::latex2html(outreach.date) << ". ";// << endl;
      oss << "<br>" << endl;
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
      if(outreach._isinvited && outreach.type=="PRESENTATION_TALK") oss << "{\\bf Invited talk:}";
      if(outreach._isinvited && outreach.type=="PRESENTATION_SEMINAR") oss << "{\\bf Invited seminar:}";
      if(outreach._isinvited && outreach.type=="PRESENTATION_COLLOQUIUM") oss << "{\\bf Invited colloquium:}";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PLENARY") oss << "{\\bf Plenary Speaker:}";
      if(outreach._isinvited && outreach.type=="PRESENTATION_KEYNOTE") oss << "{\\bf Keynote Speaker:}";
      if(outreach._isinvited && outreach.type=="PRESENTATION_TUTORIAL") oss << "{\\bf Tutorial:}";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PANELIST") oss << "{\\bf Invited panelist:}";
      oss << endl;
      // 3rd line AUTHORS
      if(0) {
	for(uint iauth=0;iauth<outreach.vauthor.size();iauth++)
	  oss << outreach.vauthor.at(iauth) << ", ";
	oss << endl;
      }
      oss << "{\\it " << outreach.title << "}; " << endl;
      oss << "" << outreach.place << ", " << endl;
      oss << "" << outreach.date << ". ";// << endl;
    }
  }
  return oss;
}


// ******************************************************************************************************************************************************
void voutreach_print_publications_EXTRA(ostream& oss);

vector<_outreach> voutreach_global_list;
uint voutreach_global_max_year=0,voutreach_global_min_year=9999;


// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// ARTICLES
// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

uint voutreach_remove_duplicate(vector<_outreach>& voutreach) {
  // cerr << voutreach.size() << endl;
  for(uint i=0;i<voutreach.size();i++)
    for(uint j=i+1;j<voutreach.size();j++)
      if(i!=j && voutreach.at(j).year==voutreach.at(i).year) { // same year
	// cerr << "found same year: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
	if(voutreach.at(j).vauthor.size()==voutreach.at(i).vauthor.size()) { // same number of authors
	  // cerr << "found same vauthor.size: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
	  if(voutreach.at(j).vauthor.at(0)==voutreach.at(i).vauthor.at(0)) { // same first author
	    // cerr << "found same vauthor.at(0): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
	    if(voutreach.at(j).vauthor.at(voutreach.at(j).vauthor.size()-1)==voutreach.at(i).vauthor.at(voutreach.at(i).vauthor.size()-1)) { // same last author
	      // cerr << "found same vauthor.at(N): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
	      vector<string> tokensi;aurostd::string2tokens(voutreach.at(i).title,tokensi," ");
	      vector<string> tokensj;aurostd::string2tokens(voutreach.at(j).title,tokensj," ");
	      if(tokensj.size()==tokensi.size()) { // same year
		// cerr << "found same title.size: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
		// cerr << "removing (i,j)=(" << i << "," << j << ")  title1=" << voutreach.at(j).title << " | " << voutreach.at(i).title << endl;
		if(i!=j) {voutreach.erase(voutreach.begin()+j);i=0;j=0;}
		//	cerr << "Article: found duplicate (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
	      }
	    }
	  }
	}
      }
  return voutreach.size();
}


// ******************************************************************************************************************************************************
void voutreach_print_publications_EXTRA(ostream& oss) {
  oss << "<! *************************************************************************************************************>" << endl;
  oss << "<! *************************************************************************************************************>" << endl;
  oss << "" << endl;
  oss << "<img border=\"0\" width=100% height=3 src=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/auro/images/line.gif>" << endl;
  oss << "" << endl;
  // **********************************************************************************************
  if(1) { // PATENTS
    oss << "<li><span class=\"pubYears\">Patents</span></li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubAuthors\">G. Ceder, C. Fischer, K. Tibbetts, D. Morgan, and S. Curtarolo</span>, " << endl;
    oss << "<span class=\"pubTitle\">Systems and Methods for predicting materials properties</span>,"  << endl;
    oss << "<span class=\"pubJournal\">US Patent #7292958 (2007)</span>. [<a HREF=" << WEB_PDF << "PAT7292958.pdf>pdf</a>]" << endl;
    oss << "</li>" << endl;
  }
  // **********************************************************************************************
  if(1) { // THESES
    oss << "<li><span class=\"pubYears\">Theses</span></li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubNumber\">4. </span><span class=\"pubAuthors\">S. Curtarolo</span>," << endl;
    oss << "<span class=\"pubTitle\">Coarse-Graining and Data Mining Approaches to the Prediction of Structures and their Dynamics</span>, " << endl;
    oss << "<span class=\"pubJournal\">Ph.D., Massachusetts Institute of Technology (2003)</span>." << endl;
    oss << "[<a href=" << WEB_PDF << "scurtarolo4.pdf>pdf</a>]<br>" << endl;
    oss << "</li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubNumber\">3. </span><span class=\"pubAuthors\">S. Curtarolo</span>," << endl;
    oss << "<span class=\"pubTitle\">Adsorption Problems investigated with Computer Simulation</span>, " << endl;
    oss << "<span class=\"pubJournal\">M.S., Pennsylvania State University (1999)</span>." << endl;
    oss << "[<a href=" << WEB_PDF << "scurtarolo3.pdf>pdf</a>]<br>" << endl;
    oss << "</li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubNumber\">2. </span><span class=\"pubAuthors\">S. Curtarolo</span>," << endl;
    oss << "<span class=\"pubTitle\">Influenza della rugosita' sul prewetting di Neon su Magnesio</span>, " << endl;
    oss << "<span class=\"pubJournal\">Laurea in Fisica, Universita` di Padova (1998)</span>." << endl;
    oss << "[<a href=" << WEB_PDF << "scurtarolo2.pdf>pdf</a>] (in Italian). <br>" << endl;
    oss << "</li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubNumber\">1. </span><span class=\"pubAuthors\">S. Curtarolo</span>," << endl;
    oss << "<span class=\"pubTitle\">Approccio analitico e numerico allo studio degli adattatori dielettrici</span>, " << endl;
    oss << "<span class=\"pubJournal\">Laurea in Ingegneria Elettronica, Universita` di Padova (1995)</span>." << endl;
    oss << "[<!A href=" << WEB_PDF << "scurtarolo1.pdf>pdf, not available<!/a>] (in Italian). <br>" << endl;
    oss << "</li>" << endl;
  }
}

// ***************************************************************************
void SystemReferences(const string& system_in,vector<uint>& voutreach_wnumber) {  // ADD REFERENCES
  voutreach_wnumber.clear();
  
  string system=KBIN::VASP_PseudoPotential_CleanName(system_in);
  vector<_outreach> voutreach;
  voutreach_load(voutreach,"PUBLICATIONS");
  //   cerr << "LOADED " << voutreach.size() << " articles" << endl;
  for(uint iart=0;iart<voutreach.size();iart++)
    if(voutreach.at(iart).valloy.size()>0) 
      for(uint itoken=0;itoken<voutreach.at(iart).valloy.size();itoken++) 
	if(system==voutreach.at(iart).valloy.at(itoken)) 
	  voutreach_wnumber.push_back(voutreach.at(iart).wnumber);
}

// ***************************************************************************

void SystemReferences(const string& system_in,vector<string>& vrefs,bool AUTHOR_ETAL) {  // ADD REFERENCES
  vrefs.clear();
  vector<uint> voutreach_wnumber;
  SystemReferences(system_in,voutreach_wnumber);
  vector<_outreach> voutreach;
  voutreach_load(voutreach,"PUBLICATIONS");
  
  stringstream oss;
  
  for(uint iarticle=0;iarticle<voutreach_wnumber.size();iarticle++) {
    for(uint i=0;i<voutreach.size();i++) {
      if(voutreach_wnumber.at(iarticle)==voutreach.at(i).wnumber) {
	oss.clear();oss.str(std::string());
	// GOT THE ARTICLE
	//	if(iarticle+1<10)  oss << "<sup>&nbsp;&nbsp;" << iarticle+1 << "</sup> "; // LABEL
	//	if(iarticle+1>=10) oss << "<sup>" << iarticle+1 << "</sup> "; // LABEL
	// AUTHORS
	string authors;
	if(!AUTHOR_ETAL || voutreach.at(i).vauthor.size()<=4) { // all authors OR <=4
	  for(uint iauth=0;iauth<voutreach.at(i).vauthor.size();iauth++) {
	    authors+=voutreach.at(i).vauthor.at(iauth);
	    if(voutreach.at(i).vauthor.size()==2 && iauth==voutreach.at(i).vauthor.size()-2) authors+=" and ";
	    if(voutreach.at(i).vauthor.size()!=2 && iauth==voutreach.at(i).vauthor.size()-2) authors+=", and ";
	    if(iauth!=voutreach.at(i).vauthor.size()-2) authors+=", ";
	  }
	  //	aurostd::StringSubst(authors,"Csányi","Csanyi");
	} else { // first et al..
	  authors+=voutreach.at(i).vauthor.at(0)+" et al., ";
	}
	oss << authors;// << endl;
	// TITLE
	string title=voutreach.at(i).title;
	title=""+title+", ";
	oss << title;
	// JOURNAL with year
	string journal=voutreach.at(i).journal;
	oss << "" << journal << ". ";
	vrefs.push_back(aurostd::html2txt(oss.str()));
	oss.clear();oss.str(std::string());
      }
    }
  }
}

// ******************************************************************************************************************************************************

bool SystemInSystems(const string& system,const string& systems) {
  vector<string> tokens;
  aurostd::string2tokens(systems,tokens,",");
  for(uint i=0;i<tokens.size();i++)
    if(system==tokens.at(i))
      return TRUE;
  return FALSE;
}

bool SystemInSystems(const string& system,const vector<string>& vsystems) {
  for(uint i=0;i<vsystems.size();i++)
    if(system==vsystems.at(i))
      return TRUE;
  return FALSE;
}

// ******************************************************************************************************************************************************

bool AlloyAlphabeticLIBRARY(const string& s) {
  // never anymore this problem
  if(s=="MoMg"||s=="NaMg"||s=="NbMg"||s=="OsMg"||s=="PbMg"||s=="RbMg"||s=="RbPd"||s=="RbPt"||s=="ReMg"||s=="RhMg"||s=="RuMg"||s=="ScMg"||s=="SiPd") return FALSE;
  if(s=="SiPt"||s=="SnMg"||s=="SnPd"||s=="SnPt"||s=="SrMg"||s=="SrPd"||s=="SrPt"||s=="TaMg"||s=="TiMg"||s=="VMg"||s=="WMg"||s=="YMg"||s=="ZnMg"||s=="ZrMg") return FALSE;
  return TRUE;
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// ARTICLES/PRESENTATIONS
// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

void voutreach_print(uint _mode,ostream& oss,string what2print) {
  aurostd::xoption vflag=XHOST.vflag_outreach;
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  uint mode=_mode;

  // ******************************************************************************************************************************************************
  // ARTICLES
  if(aurostd::substring2bool(what2print,"ARTICLE") || aurostd::substring2bool(what2print,"PUBLICATION")) {
    vector<_outreach> voutreach;
    voutreach_load(voutreach,"PUBLICATIONS");
    bool LOCAL_VERBOSE=1;//FALSE;
    cerr << "LOADED " << voutreach.size() << " articles" << endl;
    bool HTRESOURCE_MODE_PHP_PREAMBLE=FALSE;

    // if(LDEBUG) cerr << voutreach.at(0);// << endl;
    // if(LDEBUG) cerr << voutreach.at(1);// << endl;
    // if(LDEBUG) cerr << voutreach.at(2);// << endl;
    // if(LDEBUG) cerr << voutreach.at(3);// << endl;
    // if(LDEBUG) cerr << voutreach.at(4);// << endl;
    // if(LDEBUG) cerr << voutreach.at(5);// << endl;
    // if(LDEBUG) cerr << voutreach.at(6);// << endl;
    // if(LDEBUG) cerr << voutreach.at(7);// << endl;
    // if(LDEBUG) cerr << voutreach.at(8);// << endl;
    //    if(LDEBUG) exit(0); 
        
    vector<_outreach> voutreach_local;
    vector<string> vauthor,vkeyword,valloy;
    //  vector<string> vargument;
    bool flag_simple=FALSE;
  
    if(LOCAL_VERBOSE) cerr << "voutreach_print [begin]" << endl;
  
    if(LDEBUG) oss << "XHOST.vflag_control.flags" << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HTML\")=" << XHOST.vflag_control.flag("PRINT_MODE::HTML") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::TXT\")=" << XHOST.vflag_control.flag("PRINT_MODE::TXT") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::LATEX\")=" << XHOST.vflag_control.flag("PRINT_MODE::LATEX") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::YEAR\")=" << XHOST.vflag_control.flag("PRINT_MODE::YEAR") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DOI\")=" << XHOST.vflag_control.flag("PRINT_MODE::DOI") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::EXTRA\")=" << XHOST.vflag_control.flag("PRINT_MODE::EXTRA") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NUMBER\")=" << XHOST.vflag_control.flag("PRINT_MODE::NUMBER") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HYPERLINKS\")=" << XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NOTE\")=" << XHOST.vflag_control.flag("PRINT_MODE::NOTE") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NEW\")=" << XHOST.vflag_control.flag("PRINT_MODE::NEW") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DATE\")=" << XHOST.vflag_control.flag("PRINT_MODE::DATE") << endl;

   // big hack in chinatown ! if Frisco gives me nothing...
    if(XHOST.vflag_control.flag("PHP::PUBS_ALLOY")) {
      if( XHOST.vflag_control.getattachedscheme("PHP::PUBS_ALLOY")=="--print=html") {
	mode=HTRESOURCE_MODE_PHP_ALLOY;
	XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE);
	XHOST.vflag_control.flag("PRINT_MODE::TXT",FALSE);
	XHOST.vflag_control.flag("PRINT_MODE::HTML",TRUE); // DEFAULT
	valloy.push_back("AFLOW");
	valloy.push_back("AFLOWLIB");
	valloy.push_back("nmatHT");
      }
    } 

    if(mode==HTRESOURCE_MODE_NONE) {mode=HTRESOURCE_MODE_PHP_AUTHOR;} // by default
    
    if(XHOST.vflag_control.flag("CV::AUTHOR") || mode==HTRESOURCE_MODE_PHP_AUTHOR)  {
      mode=HTRESOURCE_MODE_PHP_AUTHOR;
       if(LOCAL_VERBOSE) cerr << "voutreach_print: [" << XHOST.vflag_control.flag("CV::AUTHOR") << "]" << endl;
       if(XHOST.vflag_control.flag("CV::AUTHOR")) {
	 if(LOCAL_VERBOSE) cerr << "voutreach_print: [" << XHOST.vflag_control.getattachedscheme("CV::AUTHOR") << "]" << endl;
	 aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"),vauthor,",");
	 if(vauthor.size()==0) {
	   cerr << "no authors specified" << endl;
	   exit(0);
	 }
       }
   }
 
    if(LOCAL_VERBOSE) cerr << "voutreach_print: vauthor.size()=" << vauthor.size() << endl;
 
    if(XHOST.vflag_control.flag("PHP::PUBS_ALLOY") || mode==HTRESOURCE_MODE_PHP_ALLOY)  {
      if(XHOST.vflag_control.getattachedscheme("PHP::PUBS_ALLOY")!="--print=html") {
	mode=HTRESOURCE_MODE_PHP_ALLOY;
	aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("PHP::PUBS_ALLOY"),valloy,",");
	if(valloy.size()==0) exit(0);
      }
    }
    if(XHOST.vflag_control.flag("PHP::PUBS_KEYWORD") || mode==HTRESOURCE_MODE_PHP_THRUST)  {
      mode=HTRESOURCE_MODE_PHP_THRUST;
      aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("PHP::PUBS_KEYWORD"),vkeyword,",");
      if(vkeyword.size()==0) exit(0);
      for(uint i=0;i<vkeyword.size();i++) cerr << "vkeyword.at(" << i << ")=" <<  vkeyword.at(i) << endl;//exit(0);
    }
        
    // cerr << "vauthor.size()=" <<  vauthor.size() << endl;
    // cerr << "valloy.size()=" << valloy.size() << endl;
    // cerr << "vkeyword.size()=" << vkeyword.size() << endl;
    // cerr << "mode=" << mode << endl;
    // ************************************************************************************************************************
    // ************************************************************************************************************************
    // OPERATION ON MULTIPLE VAUTHOR
    if(vauthor.size()>0) {
      //   cerr << "OPERATION ON MULTIPLE VAUTHOR" << endl;
      uint recent_articles=0;
      for(uint iauthor=0;iauthor<vauthor.size();iauthor++) {
	voutreach_local.clear();
	oss << endl;
	if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_AUTHOR) {
	  oss << "<?php" << endl;
	  if(flag_simple==FALSE)
	    oss << "if($author==\"" << vauthor.at(iauthor) << "\")" << endl;
	  oss << "{" << endl;
	  oss << "?>" << endl;
	}
	if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
	  oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	  oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	  //	oss << "% PUBLICATIONS" << endl;
	  oss << "" << endl;
	  oss << "\\section{Publications} \\label{publications}" << endl;
	  oss << "[Articles can be accessed through the embedded link. Only published, submitted and ``in press'' articles are listed. The list might be slightly out of order.]" << endl;

	  //    oss << "\\begin{list}{\\labelitemi}{\\leftmargin=1em}" << endl;
	}

	// cerr << "LOADED " << voutreach.size() << " articles" << endl;

	for(uint i=0;i<voutreach.size();i++)
	  if(aurostd::substring2bool(voutreach.at(i).vauthor,vauthor.at(iauthor)))
	    voutreach_local.push_back(voutreach.at(i));  // push authors
	
	// cerr << "LOADED_LOCAL " << voutreach_local.size() << " articles" << endl;
	voutreach_remove_duplicate(voutreach_local);
	// cerr << "LOADED_DUPLICATE " << voutreach_local.size() << " articles" << endl;
	voutreach_rsort_wnumber(voutreach_local); // sort TOP numbers first
	// cerr << "LOADED_SORTED " << voutreach_local.size() << " articles" << endl;
	
	if(mode==HTRESOURCE_MODE_PHP_AUTHOR && !XHOST.vflag_control.flag("PRINT_MODE::LATEX") && !XHOST.vflag_control.flag("PRINT_MODE::TXT")) {oss << "<br><br>" << endl;}
	if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {oss << endl;}

	// ************************************************************************************************************************
	// HTML MODE
	if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
	  //	if(mode==HTRESOURCE_MODE_PHP_AUTHOR && !XHOST.vflag_control.flag("PRINT_MODE::LATEX")&& !XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
	  if(LOCAL_VERBOSE)
	    oss << "<li><span class=\"aflow\">" << "AFLOW V" << string(AFLOW_VERSION) << " </span></li>" << endl;
	  for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
	    bool flag_year=FALSE;
	    for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
	      if(year==voutreach_local.at(i).year)
		flag_year=TRUE;
	    if(flag_year && recent_articles<AUTHOR_RECENT_ARTICLES && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
	      oss << endl;
	      if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "<li><span class=\"pubYears\">" << year << "</span></li>" << endl;
	      for(uint i=0;i<voutreach_local.size();i++) {
		string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
		if(wnumber.length()<3) wnumber="0"+wnumber;
		if(wnumber.length()<2) wnumber="0"+wnumber;
		if(voutreach_local.at(i).year==year && recent_articles<AUTHOR_RECENT_ARTICLES && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
		  // print 1 article
		  XHOST.vflag_control.flag("PRINT_MODE::HTML",TRUE);
		  oss << voutreach_local.at(i) << endl;
		  recent_articles++;
		}
	      }
	    }
	  }
	  if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_AUTHOR) {
	    oss << endl << "<?php" << endl;
	    oss << "}" << endl;
	    oss << "?>" << endl;
	  }
	  if(XHOST.vflag_control.flag("PRINT_MODE::EXTRA")) voutreach_print_publications_EXTRA(oss);
	} // HTML MODE
	// ************************************************************************************************************************

	// ************************************************************************************************************************
	// LATEX MODE
	if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
	  if(!XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\begin{itemize}" << endl;

	  for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
	    bool flag_year=FALSE;
	    for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
	      if(year==voutreach_local.at(i).year)
		flag_year=TRUE;
	    if(flag_year && year<=voutreach_global_max_year) {
	      oss << endl;
	      if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\ \\\\ { \\bf \\color{blue}{[" << year << "]}}" << endl;
	      if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\begin{itemize}" << endl;
	      for(uint i=0;i<voutreach_local.size();i++) {
		string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
		if(wnumber.length()<3) wnumber="0"+wnumber;
		if(wnumber.length()<2) wnumber="0"+wnumber;
		if(voutreach_local.at(i).year==year && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
		  // print 1 article
		  if(!XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
		    oss << "" << "\\item[$\\bullet$]{}";// << "% ART" << wnumber << " - aflow " << string(AFLOW_VERSION);
		  } else {
		    if(voutreach_local.at(i).wnumber==0) {
		      oss << "" << "\\item[{\\bf Chapter.\\,}]{}";// << "% ART" << wnumber << " - aflow " << string(AFLOW_VERSION);
		    } else {
		      oss << "" << "\\item[{\\bf "+aurostd::utype2string(voutreach_local.at(i).wnumber)+".\\,}]{}";// << "% ART" << wnumber << " - aflow " << string(AFLOW_VERSION);
		    }
		  }
		  // oss << endl;
 		  oss << voutreach_local.at(i);// << endl;
		}
	      }
	      if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\end{itemize}" << endl;
	    }
	  }
	  oss << "" << endl;
	  if(!XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\end{itemize}" << endl;
	  //    oss << "\\end{list}" << endl;
	  oss << "" << endl;
	  oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	  oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	} // LATEX MODE
	// ************************************************************************************************************************

	// ************************************************************************************************************************
	// TXT' MODE
	if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
	  for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
	    bool flag_year=FALSE;
	    for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
	      if(year==voutreach_local.at(i).year)
		flag_year=TRUE;
	    if(flag_year && year<=voutreach_global_max_year) {
	      oss << endl;
	      if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "[" << year << "]" << endl;
	      for(uint i=0;i<voutreach_local.size();i++) {
		string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
		if(wnumber.length()<3) wnumber="0"+wnumber;
		if(wnumber.length()<2) wnumber="0"+wnumber;
		if(voutreach_local.at(i).year==year && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
		  // print 1 article
		  if(XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
		    oss << "" << "["+aurostd::utype2string(voutreach_local.at(i).wnumber)+".] ";
		  }
 		  oss << voutreach_local.at(i);// << endl;
		}
	      }
	    }
	  }
	  oss << "" << endl;
	  oss << "" << endl;
	  oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	  oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	} // TXT' MODE
	// ************************************************************************************************************************
 
     } // iauthor
    }
    // ************************************************************************************************************************
    // ************************************************************************************************************************
    // OPERATION ON MULTIPLE VTHRUST
    if(vkeyword.size()>0) {
      //  cerr << "OPERATION ON MULTIPLE VTHRUST" << endl;
      uint recent_articles=0;
      voutreach_local.clear();
      for(uint ithrust=0;ithrust<vkeyword.size();ithrust++) {
	oss << endl;
	if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_THRUST) {
	  oss << "<?php" << endl;
	  if(flag_simple==FALSE)
	    oss << "if($thrust==\"" << vkeyword.at(ithrust) << "\")" << endl;
	  oss << "{" << endl;
	  oss << "?>" << endl;
	}
	for(uint i=0;i<voutreach.size();i++) {
	  //	cerr << "voutreach.at(" << i << ").vkeyword.size()=" << voutreach.at(i).vkeyword.size() << endl;
	  for(uint j=0;j<voutreach.at(i).vkeyword.size();j++)
	    if(voutreach.at(i).vkeyword.at(j)==vkeyword.at(ithrust))
	      voutreach_local.push_back(voutreach.at(i));  // push thrust
	}
      } // ithrust
      voutreach_remove_duplicate(voutreach_local);
      voutreach_sort_year(voutreach_local); // sort TOP numbers first
    
      //    if(mode==HTRESOURCE_MODE_PHP_THRUST) {oss << "<!br><!br>" << endl;}
      // ************************************************************************************************************************
      // HTML MODE
      if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
	if(mode==HTRESOURCE_MODE_PHP_THRUST) {
	  if(LOCAL_VERBOSE)
	    oss << "<li><span class=\"aflow\">" << "AFLOW V" << string(AFLOW_VERSION) << " </span></li>" << endl;
	  for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
	    bool flag_year=FALSE;
	    for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
	      if(year==voutreach_local.at(i).year)
		flag_year=TRUE;
	    if(flag_year && recent_articles<THRUST_RECENT_ARTICLES && year>=voutreach_global_max_year-THRUST_RECENT_YEARS && year<=voutreach_global_max_year) {
	      for(uint i=0;i<voutreach_local.size();i++) {
		string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
		if(wnumber.length()<3) wnumber="0"+wnumber;
		if(wnumber.length()<2) wnumber="0"+wnumber;
		if(voutreach_local.at(i).year==year && recent_articles<THRUST_RECENT_ARTICLES && year>=voutreach_global_max_year-THRUST_RECENT_YEARS && year<=voutreach_global_max_year) {
		  XHOST.vflag_control.flag("PRINT_MODE::NUMBER",FALSE);
		  oss << voutreach_local.at(i) << endl;
		  recent_articles++;
		}
	      }
	    }
	  }
	  if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_THRUST) {
	    oss << "<?php" << endl;
	    oss << "}" << endl;
	    oss << "?>" << endl;
	  }
	}
      } // HTML MODE
      // ************************************************************************************************************************
    }
    // ************************************************************************************************************************
    // ************************************************************************************************************************
    // OPERATION ON MULTIPLE VALLOY
    if(valloy.size()>0) {
      //  cerr << "OPERATION ON MULTIPLE VALLOY" << endl;
      //   uint recent_articles=0;
      voutreach_local.clear();
      for(uint ialloy=0;ialloy<valloy.size();ialloy++) {
	vector<uint> voutreach_wnumber;
	SystemReferences(KBIN::VASP_PseudoPotential_CleanName(valloy.at(ialloy)),voutreach_wnumber);
	for(uint iwnumber=0;iwnumber<voutreach_wnumber.size();iwnumber++) 
	  for(uint iarticle=0;iarticle<voutreach.size();iarticle++)
	    if(voutreach.at(iarticle).wnumber==voutreach_wnumber.at(iwnumber))
	      //	  if(voutreach.at(iart)!=65 && voutreach.at(iart)!=75)     // remove aflow and aflowlib papers
	      voutreach_local.push_back(voutreach.at(iarticle));                     // remove aflow and aflowlib papers
	// cerr << voutreach_local.size() << endl;
      
      } // ialloy
      voutreach_remove_duplicate(voutreach_local);
      voutreach_sort_year(voutreach_local); // sort TOP numbers first
    
      //    if(mode==HTRESOURCE_MODE_PHP_ALLOY) {oss << "<!br><!br>" << endl;}
      // ************************************************************************************************************************
      // PHP AND WEB MODE
      if(mode==HTRESOURCE_MODE_PHP_ALLOY) {
	if(LOCAL_VERBOSE)
	  if(XHOST.vflag_control.flag("PRINT_MODE::HTML") || mode==HTRESOURCE_MODE_PHP_AUTHOR)
	    oss << "<li><span class=\"aflow\">" << "AFLOW V" << string(AFLOW_VERSION) << " </span></li>" << endl;
	oss << "<ol class=\"reference\">" << endl;
	for(uint i=0;i<voutreach_local.size();i++) {
	  //	string wnumber=voutreach_local.at(i).wnumber;
	  //	if(wnumber.length()<2) wnumber="0"+wnumber;
	  XHOST.vflag_control.flag("PRINT_MODE::NUMBER",FALSE);
	  voutreach_local.at(i).vextra_html.clear(); // delete it
	  oss << voutreach_local.at(i) << endl;
	  // recent_articles++;
	}
	oss << "</ol>" << endl;
	if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_ALLOY) {
	  oss << "<?php" << endl;
	  oss << "}" << endl;
	  oss << "?>" << endl;
	}
      } // PHP WEB
      // ************************************************************************************************************************
    }
  
    /*
      for(uint i=0;i<voutreach_local.size();i++)
      oss << voutreach_local.at(i) << endl;
    */
  }

  // ******************************************************************************************************************************************************
  // PRESENTATIOS
  if(aurostd::substring2bool(what2print,"PRESENTATION")) {
    vector<_outreach> voutreach;
    voutreach_load(voutreach,"PRESENTATIONS");
    cerr << "LOADED " << voutreach.size() << " presentations" << endl;
    
    if(mode==HTRESOURCE_MODE_NONE) {mode=HTRESOURCE_MODE_PHP_AUTHOR;}
    
    vector<_outreach> voutreach_local;
    vector<string> vauthors;
    vauthors.push_back("Curtarolo");
    
    uint noutreach=0;
    
    // ************************************************************************************************************************
    // TXT MODE
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      for(uint iauthor=0;iauthor<vauthors.size();iauthor++) {
	voutreach_local.clear();
	
	oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	for(uint i=0;i<voutreach.size();i++)
	  if(aurostd::substring2bool(voutreach.at(i).vauthor,vauthors.at(iauthor)))
	    voutreach_local.push_back(voutreach.at(i));	
	for(uint year=MAX_YEAR_PRESENTATIONS;year>=1998;year--) {
	  bool flag_year=FALSE;
	  for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
	    if(year==voutreach_local.at(i).year) flag_year=TRUE;
	  if(flag_year) {
	    for(uint i=0;i<voutreach_local.size();i++)
	      if(voutreach_local.at(i).year==year) {
		// print 1 outreach
		oss << endl;
		oss << voutreach_local.at(i) << endl;
	      }
	  }
	}
	oss << "" << endl;
      } // iauthor
    } // TXT MODE
    // ************************************************************************************************************************
    // LATEX MODE
    if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
      for(uint iauthor=0;iauthor<vauthors.size();iauthor++) {
	voutreach_local.clear();
	
	oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	oss << "% TALKS" << endl;
	oss << "" << endl;
	oss << "\\section{Invited Talks, Seminars, Plenaries}" << endl;
	oss << "\\label{italks}" << endl;
	
	/*
	  uint number_italks=0,number_iseminars=0,number_icolloquia=0;
	  for(uint i=0;i<voutreach.size();i++) {
	  if(voutreach.at(i)._istalk) number_italks++;
	  if(voutreach.at(i)._isseminar) number_iseminars++;
	  if(voutreach.at(i)._iscolloquium) number_icolloquia++;
	  }
	  oss << "{\\bf Number of invited talks at conferences: " << number_italks << " (talks). } \\\\" << endl;
	  oss << "{\\bf Number of invited seminar and colloquia: " << number_iseminars << " (seminars), " << number_icolloquia << " (colloquia). } \\\\" << endl;
	  oss << "(the list includes the invitations already scheduled for the Fall 2011)." << endl;
	*/
	oss << "" << endl;
	oss << "\\begin{itemize}" << endl;
	// oss << "\\begin{list}{\\labelitemi}{\\leftmargin=1em}" << endl;
	
	
	for(uint i=0;i<voutreach.size();i++)
	  if(aurostd::substring2bool(voutreach.at(i).vauthor,vauthors.at(iauthor)))
	    voutreach_local.push_back(voutreach.at(i));
	
	for(uint year=MAX_YEAR_PRESENTATIONS;year>=1998;year--) {
	  bool flag_year=FALSE;
	  for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
	    if(year==voutreach_local.at(i).year) flag_year=TRUE;
	  if(flag_year) {
	    //	oss << endl << "<b><font size=\"3\"><font COLOR=blue><i>" << year << ".</i> </font></b><br><br>" << endl;
	    for(uint i=0;i<voutreach_local.size();i++)
	      if(voutreach_local.at(i).year==year) {
		// print 1 outreach
		oss << endl;
		if(vauthors.at(iauthor)=="Curtarolo") 
		  oss << "\\item[{\\bf " << voutreach_local.size()-noutreach++ << ".\\,}]{}" << endl;
		oss << voutreach_local.at(i) << endl;
	      }
	  }
	}
	oss << "\\end{itemize}" << endl;
	// 	oss << "\\end{list}" << endl;
	
	oss << "" << endl;
	oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      } // iauthor
    } // TXT MODE
    // ************************************************************************************************************************
    // HTML MODE
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      if(mode==HTRESOURCE_MODE_PHP_AUTHOR) {
	for(uint iauthor=0;iauthor<vauthors.size();iauthor++) {
	  voutreach_local.clear();
	  
	  for(uint i=0;i<voutreach.size();i++)
	    if(aurostd::substring2bool(voutreach.at(i).vauthor,vauthors.at(iauthor)))
	      voutreach_local.push_back(voutreach.at(i));
	  
	  for(uint year=MAX_YEAR_PRESENTATIONS;year>=1998;year--) {
	    bool flag_year=FALSE;
	    for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
	      if(year==voutreach_local.at(i).year) flag_year=TRUE;
	    if(flag_year) {
	      //	oss << endl << "<b><font size=\"3\"><font COLOR=blue><i>" << year << ".</i> </font></b><br><br>" << endl;
	      for(uint i=0;i<voutreach_local.size();i++)
		if(voutreach_local.at(i).year==year) {
		  // print 1 outreach
		  oss << endl;
		  if(vauthors.at(iauthor)=="Curtarolo") {
		    oss << "<font COLOR=red><b>" << voutreach_local.size()-noutreach++  << ". </b></font>";
		  }
		  XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE);
		  XHOST.vflag_control.flag("PRINT_MODE::TXT",FALSE);
		  XHOST.vflag_control.flag("PRINT_MODE::HTML",TRUE);              
		  // [OBSOLETE]		voutreach_local.at(i).print_mode="HTML";
		  oss << voutreach_local.at(i) << endl;
		}
	    }
	  }
	}
      }
    } // HTML MODE
  }  // PRESENTATIONS
}

// ******************************************************************************************************************************************************
void voutreach_print_everything(ostream& oss,const vector<string>& vitems,string msg1,string msg2,string sectionlabel) {
  bool flag_ITMZ=FALSE;
  bool flag_LIST=TRUE;
  cerr << "LOADED " << vitems.size() << " " << msg1 << endl;

  if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
    oss << "% " << msg2 << endl;
    oss << sectionlabel << endl;
    if(flag_ITMZ) oss << "\\begin{itemize}" << endl;
    if(flag_LIST) oss << "\\begin{list}{\\labelitemi}{\\leftmargin=1em}" << endl;
    for(uint i=0;i<vitems.size();i++)
      oss << "\\item{}" << " " << vitems.at(i) << "  % " << msg2 << endl;
    if(flag_ITMZ) oss << "\\end{itemize}" << endl;
    if(flag_LIST) oss << "\\end{list}" << endl;
    //    oss << "" << endl;
  }
}

// ******************************************************************************************************************************************************
vector<_outreach> voutreach_presentations;
vector<_outreach> voutreach_publications;
int voutreach_call=0;
uint voutreach_load(vector<_outreach>& voutreach,string what2print) {
  //  cerr << voutreach_call++ << endl;
  // check cache
  voutreach.clear();
  vector<string> vlabel;
  if(aurostd::substring2bool(what2print,"ARTICLE") || aurostd::substring2bool(what2print,"PUBLICATION")) {
    if(voutreach_publications.size()) {
      for(uint i=0;i<voutreach_publications.size();i++)
	voutreach.push_back(voutreach_publications.at(i));
      return voutreach.size();
    }  
  }
  if(aurostd::substring2bool(what2print,"PRESENTATION")) {
    if(voutreach_presentations.size()) {
      for(uint i=0;i<voutreach_presentations.size();i++)
	voutreach.push_back(voutreach_presentations.at(i));
      return voutreach.size();
    }  
  }
  // no cache
  // <b><font size="3"><font COLOR=blue><i>2010-current.</i> </font></b><br><br>
  _outreach ptmp;
  
  string cv2open="f144468a7ccc2d3a72ba44000715efdb";
  //  cerr << XHOST.vflag_control.getattachedscheme("CV::AUTHOR") << endl;
  if(aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="CURTAROLO" || aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="SCURTAROLO") cv2open="f144468a7ccc2d3a72ba44000715efdb";
  if(aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="OSES" || aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="COSES") cv2open="d0f1b0e47f178ae627a388d3bf65d2d2";
  if(aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="TOHER" || aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="CTOHER") cv2open="decf00ca3ad2fe494eea8e543e929068";

  vector<string> vpres,ktokens,tokens;
  if(aurostd::substring2bool(what2print,"ARTICLE") || aurostd::substring2bool(what2print,"PUBLICATION")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"PRESENTATION")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"EDUCATION")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"RESEARCH")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"ACADEMIC")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"SERVICEOUTSIDE")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"SERVICEINSIDE")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"TEACHING")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"ADVISING")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"PATENTS")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"PRESS")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"AWARDS")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  
  string iline,jline,kline;
  uint i,j,k,l;
  
  for(i=0;i<vpres.size();i++) {
    //  cout << vpres.at(i) << endl;
 
    iline=vpres.at(i);
    iline=aurostd::StringSubst(iline," ","");
    if(iline=="OBJECT={") { // found an object
      for(j=i;i<vpres.size();j++) {
	jline=vpres.at(j);
	jline=aurostd::StringSubst(jline," ","");
	if(jline.at(0)=='}') break;
      }
      // now I have an object between i+1 and j-1
      ptmp.clear();
      for(k=i+1;k<j;k++) {
	kline=vpres.at(k);
	kline=aurostd::StringSubst(kline," ","");
	aurostd::string2tokens(kline,ktokens,"=");
	aurostd::string2tokens(vpres.at(k),tokens,"=");
	//	cerr << tokens.at(1) << endl;
	string object;
	if(tokens.size()>0 && ktokens.size()>0) {
	  object=tokens.at(1);
	  for(l=2;l<tokens.size();l++)
	    object+="="+tokens.at(l);
	  //	  cerr << object << endl;
	  for(l=0;l<5&&object.size();l++) { // cleanup
	    if(object.at(0)==' ' || object.at(0)=='\"')
	      object=object.substr(1,object.length());
	    if(object.at(object.length()-1)==' ' || object.at(object.length()-1)==';' || 
	       object.at(object.length()-1)=='"') object=object.substr(0,object.length()-1);
	  }
	  //cerr << object << endl;
	  if(ktokens.at(0)=="YEAR" || ktokens.at(0)=="year") ptmp.year=aurostd::string2utype<uint>(ktokens.at(1)); // check year
	  if(ktokens.at(0)=="TYPE" || ktokens.at(0)=="type") { // check type
	    if(object=="TALK" || object=="PRESENTATION_TALK") {ptmp.type="PRESENTATION_TALK";ptmp._isinvited=TRUE;}
	    if(object=="SEMINAR" || object=="PRESENTATION_SEMINAR") {ptmp.type="PRESENTATION_SEMINAR";ptmp._isinvited=TRUE;}
	    if(object=="COLLOQUIUM" || object=="PRESENTATION_COLLOQUIUM") {ptmp.type="PRESENTATION_COLLOQUIUM";ptmp._isinvited=TRUE;}
	    if(object=="KEYNOTE" || object=="PRESENTATION_KEYNOTE") {ptmp.type="PRESENTATION_KEYNOTE";ptmp._isinvited=TRUE;}
	    if(object=="PLENARY" || object=="PRESENTATION_PLENARY") {ptmp.type="PRESENTATION_PLENARY";ptmp._isinvited=TRUE;}
	    if(object=="TUTORIAL" || object=="PRESENTATION_TUTORIAL") {ptmp.type="PRESENTATION_TUTORIAL";ptmp._isinvited=TRUE;}
	    if(object=="PANELIST" || object=="PRESENTATION_PANELIST") {ptmp.type="PRESENTATION_PANELIST";ptmp._isinvited=TRUE;}
	    if(object=="CONTRIBUTED" || object=="PRESENTATION_CONTRIBUTED") ptmp.type="PRESENTATION_CONTRIBUTED";
	    if(object=="POSTER" || object=="PRESENTATION_POSTER") ptmp.type="PRESENTATION_POSTER";
	    if(object=="ARTICLE" || object=="article") ptmp.type="ARTICLE";
	    if(object=="LABEL" || object=="label") ptmp.type="LABEL";
	    if(object=="PUBLICATION" || object=="publication") ptmp.type="ARTICLE";
	    if(object=="EDUCATION" || object=="education") ptmp.type="EDUCATION";
	    if(object=="RESEARCH" || object=="research") ptmp.type="RESEARCH";
	    if(object=="ACADEMIC" || object=="academic") ptmp.type="ACADEMIC";
	    if(object=="SERVICEOUTSIDE" || object=="serviceoutside") ptmp.type="SERVICEOUTSIDE";
	    if(object=="SERVICEINSIDE" || object=="serviceinside") ptmp.type="SERVICEINSIDE";
	    if(object=="TEACHING" || object=="teaching") ptmp.type="TEACHING";
	    if(object=="ADVISING" || object=="advising") ptmp.type="ADVISING";
	    if(object=="PATENTS" || object=="patents") ptmp.type="PATENTS";
	    if(object=="PRESS" || object=="press") ptmp.type="PRESS";
	    if(object=="AWARDS" || object=="awards") ptmp.type="AWARDS";
	  }
	  if(ktokens.at(0)=="TITLE" || ktokens.at(0)=="title") ptmp.title=object; // check title
	  if(ktokens.at(0)=="JOURNAL" || ktokens.at(0)=="journal") ptmp.journal=object; // check journal
	  if(ktokens.at(0)=="LINK" || ktokens.at(0)=="link") ptmp.link=object; // check link
	  if(ktokens.at(0)=="ARXIV" || ktokens.at(0)=="arxiv") ptmp.arxiv=object; // check arxiv
	  if(ktokens.at(0)=="SUPPLEMENTARY" || ktokens.at(0)=="supplementary") ptmp.supplementary=object; // check supplementary
	  if(ktokens.at(0)=="BIBTEX" || ktokens.at(0)=="bibtex") ptmp.bibtex=object; // check bibtex
	  if(ktokens.at(0)=="PLACE" || ktokens.at(0)=="place") ptmp.place=object; // check place
	  if(ktokens.at(0)=="DATE" || ktokens.at(0)=="date") ptmp.date=object; // check date
	  if(ktokens.at(0)=="HOST" || ktokens.at(0)=="host") ptmp.host=object; // check host
	  if(ktokens.at(0)=="ABSTRACT" || ktokens.at(0)=="abstract") ptmp.abstract=object; // check abstract
	  if(ktokens.at(0)=="PDF" || ktokens.at(0)=="pdf") ptmp.pdf=object; // check pdf
	  if(ktokens.at(0)=="DOI" || ktokens.at(0)=="doi") ptmp.doi=object; // check doi
	  if(ktokens.at(0)=="NEW" || ktokens.at(0)=="new") ptmp.newflag=TRUE; // check doi
	  if(ktokens.at(0)=="NEWFLAG" || ktokens.at(0)=="newflag") ptmp.newflag=TRUE; // check doi
	  if(ktokens.at(0)=="WNUMBER" || ktokens.at(0)=="wnumber") ptmp.wnumber=aurostd::string2utype<uint>(object); // check wnumber  
	  if(ktokens.at(0)=="AUTHOR" || ktokens.at(0)=="author") aurostd::string2tokensAdd(object,ptmp.vauthor,",");
	  if(ktokens.at(0)=="EXTRA_HTML" || ktokens.at(0)=="extra_html") ptmp.vextra_html.push_back(object);
	  if(ktokens.at(0)=="EXTRA_LATEX" || ktokens.at(0)=="extra_latex") ptmp.vextra_latex.push_back(object);
	  if(ktokens.at(0)=="KEYWORD" || ktokens.at(0)=="keyword") aurostd::string2tokensAdd(object,ptmp.vkeyword,",");
	  if(ktokens.at(0)=="SPONSOR" || ktokens.at(0)=="sponsor") aurostd::string2tokensAdd(object,ptmp.vsponsor,",");
	  if(ktokens.at(0)=="ALLOY" || ktokens.at(0)=="alloy") aurostd::string2tokensAdd(object,ptmp.valloy,",");
	  if(ktokens.at(0)=="LABEL" || ktokens.at(0)=="label") aurostd::string2tokensAdd(object,vlabel,",");
	}
      }
      
      if(ptmp.type=="") {
	for(k=i;k<=j;k++) 
	  cerr << "entry=[" << vpres.at(k) << "]" << endl;
      } else {
	if((aurostd::substring2bool(ptmp.type,"ARTICLE") || aurostd::substring2bool(ptmp.type,"PUBLICATION")) && (aurostd::substring2bool(what2print,"ARTICLE") || aurostd::substring2bool(what2print,"PUBLICATION")))
	  voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"PRESENTATION") && aurostd::substring2bool(what2print,"PRESENTATION")) {
	  if(ptmp.vauthor.size()==0) ptmp.vauthor.push_back("S. Curtarolo");  // for stefano CV
	  voutreach.push_back(ptmp);
	}
 	if(aurostd::substring2bool(ptmp.type,"EDUCATION") && aurostd::substring2bool(what2print,"EDUCATION")) voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"RESEARCH") && aurostd::substring2bool(what2print,"RESEARCH")) voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"ACADEMIC") && aurostd::substring2bool(what2print,"ACADEMIC")) voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"SERVICEOUTSIDE") && aurostd::substring2bool(what2print,"SERVICEOUTSIDE")) voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"SERVICEINSIDE") && aurostd::substring2bool(what2print,"SERVICEINSIDE")) voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"TEACHING") && aurostd::substring2bool(what2print,"TEACHING")) voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"ADVISING") && aurostd::substring2bool(what2print,"ADVISING")) voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"PATENTS") && aurostd::substring2bool(what2print,"PATENTS")) voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"PRESS") && aurostd::substring2bool(what2print,"PRESS")) voutreach.push_back(ptmp);
	if(aurostd::substring2bool(ptmp.type,"AWARDS") && aurostd::substring2bool(what2print,"AWARDS")) voutreach.push_back(ptmp);
      }
    }
  }
  // cerr << "vlabel.size()=" << vlabel.size() << endl; for(uint i=0;i<vlabel.size();i++)   cerr << vlabel.at(i) << endl; exit(0);
  
  for(uint i=0;i<voutreach.size();i++) fixlabel(vlabel,voutreach.at(i).vauthor);
  for(uint i=0;i<voutreach.size();i++) fixlabel(vlabel,voutreach.at(i).vextra_html);
  for(uint i=0;i<voutreach.size();i++) fixlabel(vlabel,voutreach.at(i).vextra_latex);
  if(aurostd::substring2bool(ptmp.type,"ARTICLE") || aurostd::substring2bool(ptmp.type,"PUBLICATION"))  
    voutreach_remove_duplicate(voutreach);  // talks can be duplicated
  // SAVE the global one
  if(voutreach_global_list.size()==0) voutreach_global_list=voutreach;
  // CHECK MAX YEAR
  for(uint i=0;i<voutreach.size();i++) if(voutreach_global_max_year<voutreach.at(i).year) voutreach_global_max_year=voutreach.at(i).year;
  // CHECK MIN YEAR
  for(uint i=0;i<voutreach.size();i++) if(voutreach_global_min_year>voutreach.at(i).year) voutreach_global_min_year=voutreach.at(i).year;
  // DONE

  // save cache
  if(aurostd::substring2bool(ptmp.type,"ARTICLE") || aurostd::substring2bool(ptmp.type,"PUBLICATION")) {
    voutreach_publications.clear();
    for(uint i=0;i<voutreach.size();i++) {
      voutreach_publications.push_back(voutreach.at(i));
      //     cerr << voutreach.at(i) << endl;
    }
  } 
  if(aurostd::substring2bool(ptmp.type,"PRESENTATION")) {
    voutreach_presentations.clear();
    for(uint i=0;i<voutreach.size();i++)
      voutreach_presentations.push_back(voutreach.at(i));
  }
  //  exit(0);
  return voutreach.size();
}


// ******************************************************************************************************************************************************
void HT_CHECK_GRANTS(ostream& oss) {//,const vector<string>& vitems,string msg1,string msg2,string sectionlabel) {
  aurostd::xoption vflag=XHOST.vflag_outreach;
  vector<_outreach> voutreach;
  voutreach_load(voutreach,"PUBLICATIONS");
  if(!vflag.flag("GRANTS")) {cerr << "no grants" << endl;exit(0);}
  cerr << "LOADED " << voutreach.size() << " " << endl;
  string grant=vflag.getattachedscheme("GRANTS");
  //  if(mode==HTRESOURCE_MODE_NONE) {mode=HTRESOURCE_MODE_PHP_AUTHOR;} // by default

  if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
    for(uint i=0;i<voutreach.size();i++) {
      bool found=FALSE;
      for(uint j=0;j<voutreach.at(i).vsponsor.size()&&found==FALSE;j++)
	if(aurostd::substring2bool(voutreach.at(i).vsponsor.at(j),grant)) found=TRUE;
      if(found) {
	XHOST.vflag_control.flag("PRINT_MODE::TXT",TRUE);
	// [OBSOLETE]	voutreach.at(i).print_mode="TXT";
	oss << voutreach.at(i) << endl;// << endl;
      }
    } 
  }
}

// ******************************************************************************************************************************************************
bool ProcessPhpLatexCv(void) {
  aurostd::xoption vflag=XHOST.vflag_outreach;

  bool Arun=FALSE;
  vector<_outreach> voutreach;
  if(!XHOST.vflag_control.flag("PRINT_MODE::HTML") && 
     !XHOST.vflag_control.flag("PRINT_MODE::LATEX") && 
     !XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
    XHOST.vflag_control.flag("PRINT_MODE::HTML",TRUE);
    //    cerr << "ERROR ProcessPhpLatexCv: a print mode must be specified: --print=latex|txt|html" << endl;
    //  exit(0);
  }
  if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE);XHOST.vflag_control.flag("PRINT_MODE::TXT",FALSE);}
  if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {XHOST.vflag_control.flag("PRINT_MODE::TXT",FALSE);XHOST.vflag_control.flag("PRINT_MODE::HTML",FALSE);}
  if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {XHOST.vflag_control.flag("PRINT_MODE::HTML",FALSE);XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE);}

  if(!Arun && XHOST.vflag_control.flag("PHP::CENTER_MISSION")) {
    center_print(HTRESOURCE_MODE_PHP_AUTHOR,cout);
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("PHP::PUBS_ALLOY")) {
    voutreach_print(HTRESOURCE_MODE_PHP_ALLOY,cout,"ARTICLES");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("PHP::PUBS_KEYWORD")) {
    voutreach_print(HTRESOURCE_MODE_PHP_THRUST,cout,"ARTICLES");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("PHP::ITALKS")) {
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR,cout,"PRESENTATIONS");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::PUBS")) {
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR,cout,"ARTICLES");
    Arun=TRUE;
  }
  if(!Arun && vflag.flag("GRANTS")) {
    HT_CHECK_GRANTS(cout);
    Arun=TRUE;
  } 
  if(!Arun && XHOST.vflag_control.flag("CV::ITALKS")) {
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR,cout,"PRESENTATIONS");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::ACADEMIC")) {
    voutreach_load(voutreach,"ACADEMIC");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Academic","ACADEMIC","\\section{Academic Positions}\\label{academic_positions}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::RESEARCH")) {
    voutreach_load(voutreach,"RESEARCH");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Research","RESEARCH","\\section{Research Experience}\\label{research_experience}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::EDUCATION")) {
    voutreach_load(voutreach,"EDUCATION");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Education","EDUCATION","\\section{Education}\\label{education}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::TEACHING")) {
    voutreach_load(voutreach,"TEACHING");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Teaching","TEACHING","\\section{Teaching Experience}\\label{teaching_experience}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::ADVISING")) {
    voutreach_load(voutreach,"ADVISING");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Advising","ADVISING","\\section{Advising Experience}\\label{advising_experience}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::AWARDS")) {
    voutreach_load(voutreach,"AWARDS");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Awards","AWARDS","\\section{Awards and Honors}\\label{awards}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::PRESS")) {
    voutreach_load(voutreach,"PRESS");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Press","PRESS","\\section{Press and news releases}\\label{pressreleases}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::PATENTS")) {
    voutreach_load(voutreach,"PATENTS");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Patents","PATENTS","\\section{Patents}\\label{patents}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::SERVICE_OUTSIDE")) {
    voutreach_load(voutreach,"SERVICEOUTSIDE");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"ServiceOutside","SERVICE OUTSIDE","\\section{Outreach and Professional Activities}\\label{outreach}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::SERVICE_INSIDE")) {
    voutreach_load(voutreach,"SERVICEINSIDE");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"ServiceInside","SERVICE INSIDE","\\section{Duke University - Academic Service Activities}\\label{duke_service}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::AUTHOR")) { // something need to be specified
    // [OBSOLETE]   XHOST.vflag_control.flag("PRINT_MODE::HTML",TRUE);  // override
    // [OBSOLETE]    XHOST.vflag_control.flag("PRINT_MODE::HTML",FALSE); // override
    // [OBSOLETE]    XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE); // override
    // [OBSOLETE]  print_mode="HTML"; // override
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR,cout,"ARTICLES");
    Arun=TRUE;
  }

  if(Arun) exit(0);
  return TRUE;
}


// ******************************************************************************************************************************************************

uint voutreach_sort_year(vector<_outreach>& voutreach) {
  if(voutreach.size()>0)
    if(voutreach.at(0).year!=0) 
      sort(voutreach.begin(),voutreach.end(),_sort_outreach_outreach_year_());
  return voutreach.size();
}
uint voutreach_rsort_year(vector<_outreach>& voutreach) {
  if(voutreach.size()>0)
    if(voutreach.at(0).year!=0) 
      sort(voutreach.begin(),voutreach.end(),_rsort_outreach_outreach_year_());
  return voutreach.size();
}
uint voutreach_sort_wnumber(vector<_outreach>& voutreach) {
  if(voutreach.size()>0)
    if(voutreach.at(0).wnumber!=0) 
      sort(voutreach.begin(),voutreach.end(),_sort_outreach_outreach_wnumber_());
  return voutreach.size();
}
uint voutreach_rsort_wnumber(vector<_outreach>& voutreach) {
  if(voutreach.size()>0) 
    if(voutreach.at(0).wnumber!=0) 
      sort(voutreach.begin(),voutreach.end(),_rsort_outreach_outreach_wnumber_());
  return voutreach.size();
}

// ******************************************************************************************************************************************************
void center_print(uint mode, ostream& oss) {
  aurostd::xoption vflag=XHOST.vflag_outreach;

  if(mode) {;} // dummy load
  
  if(XHOST.vflag_control.flag("PHP::CENTER_MISSION")) {
    //  if(mode=print_mode && mode=HTRESOURCE_MODE_LATEX)
    oss << "  <br>" << endl;
   }
}

#endif //   _AFLOW_WEB_INTERFACE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
