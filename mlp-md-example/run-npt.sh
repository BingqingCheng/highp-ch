for prefix in ch4-2.00g c-1 h-1  C1H1 C1H2 C1H3 C2H1 C3H1; do

for t in 1000 2000 4000 6000 8000; do

for p in 50 100 200 400 600 800; do 

	if [ ! -e ${prefix}-P-$p-T-$t.lmp ]; then
	echo $prefix $p $t
		sed -e "s/TEMPERATURE/$t/" -e "s/PRESSURE/$p/" -e "s/PREFIX/$prefix/" npt.lmp > ${prefix}-P-$p-T-$t.lmp
	sed -e "s/TEMPERATURE/$t/g" -e "s/PRESSURE/$p/g" -e "s/PREFIX/${prefix}/g" single-npt.sh | sbatch
	fi
done
done

done
