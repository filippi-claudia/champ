#!/bin/zsh
rm -rf ~/Dropbox/Cartesius/box.tar
tar --exclude='compile.zsh' -cvf box.tar ./
ssh  plopez@cartesius.surfsara.nl 'rm -rf /home/plopez/Programs/champ/*'
scp -r box.tar plopez@cartesius.surfsara.nl:/home/plopez/Programs/champ
sleep 3 
echo
echo "Copy done!"
echo
ssh plopez@cartesius.surfsara.nl  'cd /home/plopez/Programs/champ/ ; tar -xvf box.tar ;  \
                                   rm -rf box.tar ; \
                                   cp ../compile_121919.sh  ./compile.sh ; ./compile.sh'
rm -rf box.tar