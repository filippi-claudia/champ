#!/usr/bin/gawk -f 
# common block check 
# assumes output from gnu-nm 
#  nm -g --size-sort -t d --print-file-name  
#$Revision: 1.1 $
BEGIN{
    print "** gmemchk **"; 
# size for 'bigblocks'
  if(bigsz=="") bigsz=1024*1024;
# alignment size 
  if(alsz=="") alsz=8;
  err=0;
}
# NF != 3 {error("format error in nm output"); }
# ignore lines with wrong number of fields
NF == 3 { 
  if((nfld=split($1,wtmp,":"))<2){ 
    error("format error in nm output");
  } else {
    fle=wtmp[1];
    sz=wtmp[nfld];
    ty=$2;
    sym=$3;
    if(tolower(ty)=="c"){
      c=++count[sym];
      files[sym,c]=fle;
      sizes[sym,c]=sz;
    }
  }
}

END{
  if(err) exit(1);
  tmem=0;
  bigmem=0;
  bcount=0;
  ticon=0;
  tal=0;
  for( s in count ){
    sz=sizes[s,1];
    tmem+=sz;
    if(sz>bigsz){
      bigmem+=sz;
      bigblocks[s]=sz;
    }
    if(align){
	if((sz%alsz)!=0){
	    tal++;
	    printf "** possible alignment problem for symbol %s\n",s;
	    printf "   file=%s\n",files[s,1];
	    printf "   size=%d  align=%d mod=%d\n",sz,alsz,sz%alsz;
	}
    }
    incons=0;
    for(i=2;i<=count[s];i++){
      if(sizes[s,i]!=sz){
	incons=1;
      }
    }
    if(incons==1){
      ticon++;
      printf "\n ** inconsistent sizes : symbol=%s \n",s;
      for(i=1;i<=count[s];i++){
	printf "size=%d in file %s \n",sizes[s,i],files[s,i];
      }
    }
  }
  meg=(1024*1024);
  printf "\n total size is %d (%d M)\n",tmem,tmem/meg;
  printf " blocks > %d (%d M) hold %d (%d M)\n",bigsz,bigsz/meg,bigmem,bigmem/meg;
  printf " -- big blocks --\n";
  printf "symbol \t size \t ( Meg )\n"
  for( i in bigblocks ){
    printf "%s \t %d \t %d \n",i,bigblocks[i],bigblocks[i]/meg;
  }
  if(ticon){
      printf "\n summary: inconsistencies %d \n",ticon;
  }else{
      print "no inconsistencies found";
  }
  if(align){ 
      printf " alignment problems %d\n",tal;
  } else {
      printf "no check for alignment problems\n";
      printf "use -v align=1 \n";
  } 
}

function error(msg){
  printf "*** ERROR ***\n" > "/dev/stderr";
  printf "file= %s line= %d \n",FILENAME,FNR  > "/dev/stderr";
  printf "%s\n",msg > "/dev/stderr";
  err=1;
  exit(1);
}
