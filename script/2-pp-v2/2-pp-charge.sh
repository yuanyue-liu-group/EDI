Ngrid=180
pwd=`pwd`
#################################
#2-charge
#################################
#for i in 3DLFAns    3DLFAs     3DLFAes  3DcutLFAes 3DcutLFAns 3DcutLFAs  2DLFAes    2DLFAns    2DLFAs     
for i in      3DLFAes  3DcutLFAes 
do
rm -r  mcharge$i
mkdir mcharge$i
mkdir mcharge$i/out
cp plot-charge.sh mu-charge.py mcharge$i
ln -s $pwd/kq$Ngrid.dat mcharge$i
ln -s $pwd/v$Ngrid.dat mcharge$i


for k in out/k[1-9]*.ml;
do
j=`basename $k |cut -d '.' -f 1`
ln -s $pwd/out/$j.mcharge$i $pwd/mcharge$i/out/$j.mcharge
done

cd  mcharge$i
python mu-charge.py
bash<plot-charge.sh
cd $pwd

echo $i 'plot done'
done
#################################




