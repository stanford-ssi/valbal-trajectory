#run this script from top level dir
if [ -d ignored/balloons-VALBAL ]; then
 	cd ignored/balloons-VALBAL
	git pull
else
	git clone https://github.com/stanford-ssi/balloons-VALBAL ignored/balloons-VALBAL
fi