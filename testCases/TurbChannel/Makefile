run:
	mpirun -np 8 ../../incompact3d Turbchannel.in

vis:
	~/incompact3d/visualize Turbchannel.pprc

clean:
	rm -rf *.dat ux* uy* uz* pp* vort* qcrit* Ftx* Fty* Ftz* nuSGS* probe* ALM *.xdmf
