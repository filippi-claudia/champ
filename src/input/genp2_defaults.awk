#!/usr/bin/awk -f 
#$Revision: 1.2 $
# from genp2 Revision 1.5
# additional support for default args 
# makes srs p2init and p2call from list of commands/subroutines
# default output-file p2prog.F
# input: 
#    command [sr|*|end] arglist
# or # comment
# arglist: [i | d | a | inp | defset | <lit> ] ... 
# init values for defaults : type=value (without blanks) 
# i=interger,d=double,a=character*64,inp=input-unit
# defset : current set 
# lit = ['string'| int-number | double-number]
BEGIN{
# defaults
  seterrmesg();
# dimensions (->p2_dim.inc)
  MXDEF=128;
  MXIDL=64;
  print "** genp2 (with defaults)  $Revision: 1.2 $ **"
  if(ofile=="") ofile="p2prog.F";
  indent="      ";
  splt="     $    ";
  scol=70;
  cc=0;
  lc=1;
  errorflag=0;
  endcmd="end";
# dstruc: elseif or goto
#  dstruc="elseif";
  dstruc="goto";
}
$0 ~ /^#/{next;}
$0 ~ /^[ ]*$/{next;}
$0 != ""{
  if(NF<2) {error(1)}
  if(! ($1 in dnames)){
    cc++;
    name[cc]=$1;
    dnames[$1]=1;
  } else {
    printf "name=%s\n",$1 > "/dev/stderr";
    error(9);
  }
  if($2=="*"){sname[cc]=name[cc]}else{sname[cc]=$2}
  fmt[cc]="";
  ncarg=0;
  if(NF>2){
    nargs[cc]=NF-2;
    ndefval[cc]=0;
    firstdef[cc]=0;
    firstarg[cc]=0;
    for(i=3;i<=NF;i++)
      {
	atmp=$i; 
	if((ix=index(atmp,"="))){
	  argstr=substr(atmp,1,ix-1);
          dstr=substr(atmp,ix+1);
          de=1;
	} else {
          argstr=atmp;
          de=0;
	}

	if((argstr=="i")||(argstr=="d")||(argstr=="a")){
	  if(de){
	    ndefval[cc]++;
	    defval[cc,ndefval[cc]]=dstr;
	    if(firstdef[cc]==0) { firstdef[cc]=i-2; firstarg[cc]=ncarg+1;}
	  } else {
	    if(firstdef[cc]!=0) error(4);
	  }
	} else {
           if(de) error(6);
	}
	litarg[cc,i-2]=0;
	if(argstr ~ /^'[^']*'$/){
#'
          litarg[cc,i-2]=argstr;
        } else if (argstr ~ /^[-]?[0-9]+$/){
         litarg[cc,i-2]=argstr;
        } else if (argstr ~ /^[-]?[0-9]+\.[0-9]*d[+\-]?[0-9]+$/ ){
         litarg[cc,i-2]=argstr;
        } else if(argstr=="i"){ 
          fmt[cc]=fmt[cc] "i";
          ncarg++;          
          litarg[cc,i-2]=argstr;
	}
	else if(argstr=="d"){
	  fmt[cc]=fmt[cc] "d";
           ncarg++;
           litarg[cc,i-2]=argstr;
	}
	else if(argstr=="a"){
	  fmt[cc]=fmt[cc] "a";
          litarg[cc,i-2]=argstr;
          ncarg++;
	}
       	else if(argstr=="inp"){
	  litarg[cc,i-2]=argstr;
	}
   	else if(argstr=="defset"){
	  litarg[cc,i-2]=argstr;
	}
	else {
	  print "bad arg type: ",argstr;
	  error(2);
	}
      }
  }
}

END{
  if(errorflag) exit(errorflag);
  print "C -------------------------- gen2p -------------------------------------" > ofile;
  print "C keywords from file",FILENAME >> ofile;
  print "C ----------------------------------------------------------------------" >> ofile;
  print "(",cc," lines processed)";
  p2init_h();
  p2init_b();
  p2init_t();
  p2def_h();
  p2def_b();
  p2def_t();
  p2call_h();
  p2call_b();
  p2call_t();
  print "genp2: file ",ofile," created";
}

function p2init_h(){
  pfl("subroutine p2init");
  pfl("implicit double precision (a-h,o-z)");
  pfl("include 'inc/p2_dim.inc'");
  pfl("include 'inc/p2.inc'");
  pfl("call p2ini0");
  pfl("call p2inid");
  pfl("sfile='" FILENAME "'");
}

function p2init_b(){
  pfl("nkey=" cc);
  for(k=1;k<=cc;k++){
    pfl("keys(" k ")='" name[k] "'");
    klen=length(name[k]);
    pfl("keylen(" k ")=" klen );
    na=length(fmt[k]);
    pfl("nargs(" k ")=" na);
    if(na>0){
      pfl("fmts(" k ")='" fmt[k] "'");
    } else {
      pfl("fmts(" k ")='*'");
    }
  }
}
function p2init_t(){
  pfl("end");
}


function  p2def_h(){
  pfl("subroutine p2inid");
  pfl("implicit double precision (a-h,o-z)");
  pfl("include 'inc/p2_dim.inc'");
  pfl("include 'inc/p2.inc'");
  pfl("include 'inc/p2defv.inc'");
  pfl("do i=1,MXKEY");
  pfl(" ideflt(i)=0");
  pfl(" do j=1,MXIDL");
  pfl("  idefpp(j,i)=0");
  pfl(" enddo");
  pfl("enddo");
}
function  p2def_b(){
 i0=0; f0=0; a0=0;
 has_def=0;
 for(k=1;k<=cc;k++){
   if(firstdef[k]){
     has_def=1;
     argcount=firstarg[k]-1;
     pfl("ideflt(" k ")=" firstarg[k]);
     j0=0;
     for(j=firstdef[k];j<=nargs[k];j++){
       argcount++;
       if(argcount>MXDEF) error(8);    
       if(litarg[k,j]=="i"){
        i0++; j0++;
        pfl("idefpp(" argcount "," k ")=" i0);
        pfl("idefvv(" i0 ")=" defval[k,j0]);
      } else if(litarg[k,j]=="d"){
        f0++; j0++;
        pfl("idefpp(" argcount "," k ")=" f0);
        pfl("ddefvv(" f0 ")=" defval[k,j0]);
      } else if(litarg[k,j]=="a"){
        a0++; j0++;
        pfl("idefpp(" argcount "," k ")=" a0);
        if(length(defval[k,j0])>MXIDL) error(7);
        pfl("adefvv(" a0 ")='" defval[k,j0] "'");
      } else if(litarg[k,j]=="inp"){
          error(5);
      } else if(litarg[k,j]=="defset"){
	error(10);
      } else {
        argcount--;
      }
    }
#  } else {
#     pfl("ideflt(" k ")=0");
  }
 }
}


function  p2def_t(){
 if(has_def==1){
  pfl("ip2dfl=1");
}
  pfl("end");
}


function p2call_h(){
  pfl("subroutine p2call(ikw,itmp,ftmp,is1,is2,lne,iend,MXF,iu)");
  pfl("implicit double precision (a-h,o-z)");
  pfl("include 'inc/p2etc.inc'");
  pfl("dimension itmp(MXF)");
  pfl("dimension ftmp(MXF)");
  pfl("character lne*(*)");
  pfl("dimension is1(MXF)");
  pfl("dimension is2(MXF)");
  pfl("iend=0");
}
function p2call_b(){
  ie=0;
  if(dstruc=="goto"){
    cline="goto(";
    for(k=1;k<cc;k++){
      cline=cline k ",";
    }
    cline=cline cc ") ikw";
    pfl(cline);
    pfl("call fatal('p2call: bad keyword-ID')")
  }
  for(k=1;k<=cc;k++){
    if(dstruc=="elseif"){
      if(k==1){
	pfl("if(ikw.eq.1)then");
      } else {
	pfl("elseif(ikw.eq." k ")then");
      }
    } else if( dstruc=="goto"){
      if(k<10){
	lab=" " k "   ";
      }
      else if(k<100){
	lab=" " k "  ";
      }
      else if(k<1000){
	lab=" " k " ";
      }
      else if(k<9999){
	lab=" " k ;
      } else {
	error(11);
      }
      printf("%s continue\n",lab) >> ofile;
    } else {
      error(3);
    }
    if(nargs[k]>0){
      ic=0;
      fc=0;
      ac=0;
      al="(";
      for(j=1;j<=nargs[k];j++){
        c1=litarg[k,j];      
	if(c1=="i"){
	  ic++;
	  al=al "itmp(" ic ")";
	}
	else if(c1=="d"){
	  fc++;
	  al=al "ftmp(" fc ")";
	}
	else if(c1=="a"){
	  ac++;
	  al=al "lne(is1(" ac "):is2(" ac "))";
	} 
        else if(c1=="inp"){
          al=al "iu";
        }
	else if(c1=="defset"){
          al=al "idefset";
        }
        else {
          al=al c1;
        }
	if(j<nargs[k]){
	  al=al ",";
	}else{
	  al=al ")";
	}
      }
     } else {
     	al="";
     }    
    if(sname[k]!="end"){
      pfl(" call " sname[k] al);
    } else {
      pfl(" iend=1");
      ie=1;
    } 
    if( dstruc=="goto"){
      pfl("goto 9999") ;
    }
  }
  if(ie==0){
    print "genp2: WARNING:  no end-keyword defined";
  }
}
function p2call_t(){
   if( dstruc=="elseif"){
     pfl("else");
     pfl(" stop 'p2call'");
     pfl("endif");
   } else if( dstruc=="goto"){
     printf(" 9999 continue\n") >> ofile;
   } else {
     error(3);
   }
   pfl("end ");
}
	



function pfl(txt){
  lne=indent txt;
  if(length(lne)<=scol){
    print lne >> ofile;
  }else{
    lcc=length(lne)-scol;
    ll=scol-length(splt);
    lcu=int(lcc/ll)+1;
    print substr(lne,1,scol)>>ofile;
    for(i=1;i<=lcu;i++){
      print splt substr(lne,scol+(i-1)*ll+1,ll)>>ofile;
    }
  }
}

function seterrmesg()
{
  errmesg[1]="line must have two entries at least";
  errmesg[2]="arguments must be one of i,d,a or a constant";
  errmesg[3]="dstruc must be one of elseif,goto";
  errmesg[4]="default values must continue to the last arg";
  errmesg[5]="<inp> cannot have a default value";
  errmesg[10]="<defset> cannot have a default value";
  errmesg[6]="default values only for i,d,a arguments";
  errmesg[7]="string argument too long";
  errmesg[8]="command with too many arguments";
  errmesg[9]="name multiply defined"
}
function error(X)
{
  printf "\n\n" > "/dev/stderr";
  printf "%s\n",errmesg[X] > "/dev/stderr";
  print "ERROR : " ,FILENAME ," err= ",X," at input-line ",FNR > "/dev/stderr";
  errorflag=X;
  printf "\n\n" > "/dev/stderr";
  exit(X);
}
