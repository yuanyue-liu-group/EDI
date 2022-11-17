pwd=`pwd`

#################################
#1-grep
#################################
#mkdir out
#cd ..
#for i in k[1-9]*;
#do
#grep Ml  $i/mos2.$i.out  >$pwd/out/$i.ml
#grep Mnl $i/mos2.$i.out  >$pwd/out/$i.mnl
#
#grep Mcharge3DLFAns      $i/mos2.$i.out >$pwd/out/$i.mcharge3DLFAns      
#grep Mcharge3DLFAs       $i/mos2.$i.out >$pwd/out/$i.mcharge3DLFAs      
#grep Mcharge3DcutLFAns   $i/mos2.$i.out >$pwd/out/$i.mcharge3DcutLFAns      
#grep Mcharge3DcutLFAs    $i/mos2.$i.out >$pwd/out/$i.mcharge3DcutLFAs      
#grep Mcharge2DLFAes      $i/mos2.$i.out >$pwd/out/$i.mcharge2DLFAes      
#grep Mcharge2DLFAns      $i/mos2.$i.out >$pwd/out/$i.mcharge2DLFAns     
#grep Mcharge2DLFAs       $i/mos2.$i.out >$pwd/out/$i.mcharge2DLFAs     
#
#done
#
#cd $pwd
#echo 'grep done'
#################################

#################################
#2-neutral
#################################
#mkdir mall
#python mu-all.py
#mv data.pickle muvsElimit.plt ki.plt kf-kis/ mall
#
#mkdir ml
#python mu-ml.py
#mv data.pickle muvsElimit.plt ki.plt kf-kis/ ml
#
#mkdir mnl
#python mu-mnl.py
#mv data.pickle muvsElimit.plt ki.plt kf-kis/ mnl
#echo 'python done'
#
#
#
#
#for i in ml mnl mall
#do 
#cd $i
#cp ../plot-charge.sh .
#bash <plot-charge.sh
#cd ..
#done
#echo 'plot done'
#################################


#################################
#2-charge
#################################
#for i in 3DLFAns    3DLFAs     3DcutLFAns 3DcutLFAs  2DLFAes    2DLFAns    2DLFAs     
#do
#rm -r  mcharge$i
#mkdir mcharge$i
#mkdir mcharge$i/out
#cp plot-charge.sh mu-charge.py mcharge$i
#ln -s $pwd/kq180.dat mcharge$i
#ln -s $pwd/v180.dat mcharge$i
#
#
#for k in out/k[1-9]*.ml;
#do
#j=`basename $k |cut -d '.' -f 1`
#ln -s $pwd/out/$j.mcharge$i $pwd/mcharge$i/out/$j.mcharge
#done
#
#cd  mcharge$i
#python mu-charge.py
#bash<plot-charge.sh
#cd $pwd
#
#echo 'plot done'
#done
#################################




