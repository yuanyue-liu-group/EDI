#Ngrid
ngrid=$1

# neutrla/charge type
type=$1


#################################
#2-neutral
#################################
mkdir mall
python mu-all.py
mv data.pickle muvsElimit.plt ki.plt kf-kis/ mall

mkdir ml
python mu-ml.py
mv data.pickle muvsElimit.plt ki.plt kf-kis/ ml

mkdir mnl
python mu-mnl.py
mv data.pickle muvsElimit.plt ki.plt kf-kis/ mnl
echo 'python done'




for i in ml mnl mall
do 
cd $i
cp ../plot-charge.sh .
bash <plot-charge.sh
cd ..
done
echo 'plot done'
#################################


