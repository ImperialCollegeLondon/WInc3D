run:
	mpirun -np 16 ~/WInc3D/incompact3d NREL-5MW.in

vis:
	~/WInc3D/visualize NREL-5MW.pprc

clean:
	rm -rf *.dat ux* uy* uz* pp* vort* qcrit* temp* Ftx* Fty* Ftz* nuSGS* probe* ALM *.xdmf
