run:
	mpirun -np 4 ~/incompact3d/incompact3d NREL-5MW.in

vis:
	~/incompact3d/visualize NREL-5MW.pprc

clean:
	rm -rf *.dat ux* uy* uz* pp* vort* qcrit* temp* Ftx* Fty* Ftz* nuSGS* probe* ALM *.xdmf
