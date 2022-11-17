pwd=`pwd`

#################################
#1-grep
#################################
mkdir out
cd ..
for i in k[1-9]*;
do
grep Ml  $i/mos2.$i.out  >$pwd/out/$i.ml
grep Mnl $i/mos2.$i.out  >$pwd/out/$i.mnl

grep Mcharge3DLFAns      $i/mos2.$i.out >$pwd/out/$i.mcharge3DLFAns      
grep Mcharge3DLFAs       $i/mos2.$i.out >$pwd/out/$i.mcharge3DLFAs      
grep Mcharge3DLFAes      $i/mos2.$i.out >$pwd/out/$i.mcharge3DLFAes      
grep Mcharge3DcutLFAns   $i/mos2.$i.out >$pwd/out/$i.mcharge3DcutLFAns      
grep Mcharge3DcutLFAs    $i/mos2.$i.out >$pwd/out/$i.mcharge3DcutLFAs      
grep Mcharge3DcutLFAes   $i/mos2.$i.out >$pwd/out/$i.mcharge3DcutLFAes      
grep Mcharge2DLFAns      $i/mos2.$i.out >$pwd/out/$i.mcharge2DLFAns     
grep Mcharge2DLFAs       $i/mos2.$i.out >$pwd/out/$i.mcharge2DLFAs     
grep Mcharge2DLFAes      $i/mos2.$i.out >$pwd/out/$i.mcharge2DLFAes      

done

cd $pwd
echo 'grep done'
#################################


