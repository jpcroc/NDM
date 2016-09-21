

radius=33.5
pelipsosize=20.0

width=879
height=548
antitr=0.0

grep -n -m 1 "^a(*" $1 > grep.out
#cat grep.out
nline=`sed -e "s/\:/  /g" grep.out | awk '//{x=$(0); ; printf("%i",x)}'`
echo $nline

ntot=`wc -l $1 | awk '//{x=$(0); ; printf("%i",x)}'`


nhead=`awk -v a=$nline  'BEGIN {b=a-1 ; printf("%i",b)}'`
ntail=`awk -v a=$ntot -v b=$nline  'BEGIN {c=a-b+1 ; printf("%i",c)}'`
echo $ntot $nhead $ntail

head -$nhead $1 > tmp


echo "#macro pelipso(X,Y,Z,l,R,G,B,T)" >> tmp
echo " superellipsoid{ <0.25,0.25>" >> tmp
echo "  rotate<0,18,0>" >> tmp
echo "  rotate<10,0,0>" >> tmp
echo "  scale<l,l,l>" >> tmp
echo "  translate <X,Y,Z>" >> tmp
echo "  pigment{rgbt<R,G,B,T>}" >> tmp
echo "  translucentFinish(T)" >> tmp
echo "  clip()" >> tmp
echo "  check_shadow()}" >> tmp
echo "#end" >> tmp

tail -$ntail $1 >> tmp

mv tmp $1"p"


int_radius=$(echo "$radius" | awk -F. '{print $1}')
dec_radius=$(echo "$radius" | awk -F. '{print $2}')
grep -n "\,${int_radius}\.${dec_radius}" $1"p" > grep.out
echo "We will change everthing here ...."
cat grep.out
sed -e "s/\:/  /g" grep.out > sed.out
cat sed.out | awk '//{x=$(0); ; printf("%i\n",x)}' > list.line
list=`cat list.line`
rm -f tmp.file
 for tmp in $list ; do
 echo $tmp
 sed -e "${tmp}s/a/pelipso/" $1"p" > tmp.file
 sed -e "${tmp}s/${radius}/${pelipsosize}/" tmp.file > $1"p"
 done
 
sed -e "s/rgb <0.6,0.6,0.6/rgb <0.8,0.8,0.8/g" $1"p" > tmp.file
mv tmp.file $1"p"

rm -f grep.out sed.out tmp.file list.line


sed -e s/=$1/=$1"p"/ $1".ini" > tmp
sed  "s/\(^Width=\)\(..*$\)/\1${width}/" tmp  > tmp2
sed  "s/\(^Height=\)\(..*$\)/\1${height}/" tmp2  > tmp
sed  "s/\(^Antialias_Threshold=\)\(..*$\)/\1${antitr}/" tmp  >  $1"p.ini"
