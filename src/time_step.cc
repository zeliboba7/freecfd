/************************************************************************
	
	Copyright 2007-2010 Emre Sozer

	Contact: emresozer@freecfd.com

	This file is a part of Free CFD

	Free CFD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    Free CFD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/
#include "time_step.h"

void set_time_step_options(void) {
	
	if(input.section("timemarching").get_string("integrator")=="backwardEuler") {
		time_integrator=BACKWARD_EULER;
	}

	time_step_ramp=input.section("timemarching").subsection("ramp").is_found;
	time_step_ramp_initial=input.section("timemarching").subsection("ramp").get_double("initial");
	time_step_ramp_growth=input.section("timemarching").subsection("ramp").get_double("growth");
	
	time_step_type=NONE;
	time_step_current=input.section("timemarching").get_double("stepsize");
	time_step_target=time_step_current;

	if (input.section("timemarching").get_double("stepsize").is_found) {
		time_step_type=FIXED;
	}
	if (input.section("timemarching").get_double("CFLmax").is_found) {
		if (time_step_type==FIXED) {
			cerr << "[E] Input entry timemarching -> CFLmax can't be specified together with stepsize!!" << endl;
			exit(1);
		}
		time_step_type=CFL_MAX; CFLmax=input.section("timemarching").get_double("CFLmax");
		CFLmaxTarget=CFLmax;
	}
	else if (input.section("timemarching").get_double("CFLlocal").is_found) {
		if (time_step_type==FIXED) {
			cerr << "[E] Input entry timemarching -> CFLlocal can't be specified together with stepsize!!" << endl;
			exit(1);
		}
		if (time_step_type==CFL_MAX) {
			cerr << "[E] Input entry timemarching -> CFLlocal can't be specified together with CFLmax!!" << endl;
			exit(1);
		}
		time_step_type=CFL_LOCAL; CFLlocal=input.section("timemarching").get_double("CFLlocal");
		CFLlocalTarget=CFLlocal;
	}

	// Allocate the time step variable
	dt.resize(grid.size());
	for (int gid=0;gid<grid.size();++gid) dt[gid].allocate(gid);

	return;
}

void set_pseudo_time_step_options(void) {

	// Set pseudo time step options
	ps_time_step_ramp=input.section("pseudotime").subsection("ramp").is_found;
	ps_time_step_ramp_initial=input.section("pseudotime").subsection("ramp").get_double("initial");
	ps_time_step_ramp_growth=input.section("pseudotime").subsection("ramp").get_double("growth");
	
	ps_time_step_type=NONE;
	ps_time_step_current=input.section("pseudotime").get_double("stepsize");
	ps_time_step_target=ps_time_step_current;
	
	if (input.section("pseudotime").get_double("stepsize").is_found) {
		ps_time_step_type=FIXED;
	}
	if (input.section("pseudotime").get_double("CFLmax").is_found) {
		if (ps_time_step_type==FIXED) {
			cerr << "[E] Input entry pseudotime -> CFLmax can't be specified together with stepsize!!" << endl;
			exit(1);
		}
		ps_time_step_type=CFL_MAX; ps_CFLmax=input.section("pseudotime").get_double("CFLmax");
		ps_CFLmaxTarget=ps_CFLmax;
	}
	else if (input.section("pseudotime").get_double("CFLlocal").is_found) {
		if (ps_time_step_type==FIXED) {
			cerr << "[E] Input entry pseudotime -> CFLlocal can't be specified together with stepsize!!" << endl;
			exit(1);
		}
		if (ps_time_step_type==CFL_MAX) {
			cerr << "[E] Input entry pseudotime -> CFLlocal can't be specified together with CFLmax!!" << endl;
			exit(1);
		}
		ps_time_step_type=CFL_LOCAL; ps_CFLlocal=input.section("pseudotime").get_double("CFLlocal");
		ps_CFLlocalTarget=ps_CFLlocal;
	}
	
	// Allocate the pseudo time step variable
	dtau.resize(grid.size());
	for (int gid=0;gid<grid.size();++gid) dtau[gid].allocate(gid);
		
	
	return;
}

void update_time_step_options(void) {
	input.refresh();
	input.read("timemarching");
	set_time_step_options();
	return;
}

void update_pseudo_time_step_options(void) {
	input.refresh();
	input.read("pseudotime");
	set_pseudo_time_step_options();
	return;
}

void update_time_step(int timeStep,double &time,double &max_cfl,int gid) {
	
	// TODO: What to do with ramping when restarting
	if (time_step_ramp) {
		if (timeStep==1) {
			if (time_step_type==FIXED) time_step_current=time_step_ramp_initial;
			if (time_step_type==CFL_MAX) CFLmax=time_step_ramp_initial;
			if (time_step_type==CFL_LOCAL) CFLlocal=time_step_ramp_initial;
		} else {
			if (time_step_type==FIXED) time_step_current=min(time_step_current*time_step_ramp_growth,time_step_target);
			if (time_step_type==CFL_MAX) CFLmax=min(CFLmax*time_step_ramp_growth,CFLmaxTarget);
			if (time_step_type==CFL_LOCAL) CFLlocal=min(CFLlocal*time_step_ramp_growth,CFLlocalTarget);
		}
	}
	
	double a;
	if (time_step_type==FIXED) {
		for (int gid=0;gid<grid.size();++gid) for (int c=0;c<grid[gid].cellCount;++c) dt[gid].cell(c)=time_step_current;
	} else if (time_step_type==CFL_MAX) {
		double min_dt=1.e20;
		for (int gid=0;gid<grid.size();++gid) {
			if (equations[gid]==NS) {
				time_step_current=1.e20;
				for (int c=0;c<grid[gid].cellCount;++c) {
					a=ns[gid].material.a(ns[gid].p.cell(c),ns[gid].T.cell(c));
					time_step_current=min(time_step_current,CFLmax*grid[gid].cell[c].lengthScale/(fabs(ns[gid].V.cell(c))+a));
				}
				MPI_Allreduce(&time_step_current,&time_step_current,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
				if (min_dt>time_step_current) min_dt=time_step_current;
			}
			time_step_current=min_dt;
			for (int c=0;c<grid[gid].cellCount;++c) dt[gid].cell(c)=time_step_current;
		}
	} else if (time_step_type==CFL_LOCAL) {
		// Determine time step with CFL condition
		for (int c=0;c<grid[gid].cellCount;++c) {
			a=ns[gid].material.a(ns[gid].p.cell(c),ns[gid].T.cell(c));
			dt[gid].cell(c)=CFLlocal*grid[gid].cell[c].lengthScale/(fabs(ns[gid].V.cell(c))+a);
			time_step_current=1.;
		}
	}

	max_cfl=-1.;
	double min_dt=1.e20;
	if (equations[gid]==NS) {
		for (int c=0;c<grid[gid].cellCount;++c) {
			a=ns[gid].material.a(ns[gid].p.cell(c),ns[gid].T.cell(c));
			max_cfl=max(max_cfl,(fabs(ns[gid].V.cell(c))+a)*time_step_current/grid[gid].cell[c].lengthScale);
			if (time_step_type==CFL_LOCAL) min_dt=min(min_dt,dt[gid].cell(c));
		}
		if (time_step_type==CFL_LOCAL) {
			MPI_Allreduce(&min_dt,&min_dt,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		} else {
			MPI_Allreduce(&max_cfl,&max_cfl,1, MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		}
	}
	
	if (time_step_type==CFL_LOCAL) {
		time_step_current=min_dt;
		max_cfl=CFLlocal;
	}
	time+=time_step_current;

	return;
}

void update_pseudo_time_step(int ps_step,double &max_cfl,int gid) {
	
	if (ps_time_step_ramp) {
		if (ps_step==1) {
			if (ps_time_step_type==FIXED) ps_time_step_current=ps_time_step_ramp_initial;
			if (ps_time_step_type==CFL_MAX) ps_CFLmax=ps_time_step_ramp_initial;
			if (ps_time_step_type==CFL_LOCAL) ps_CFLlocal=ps_time_step_ramp_initial;
		} else {
			if (ps_time_step_type==FIXED) ps_time_step_current=min(ps_time_step_current*ps_time_step_ramp_growth,ps_time_step_target);
			if (ps_time_step_type==CFL_MAX) ps_CFLmax=min(ps_CFLmax*ps_time_step_ramp_growth,ps_CFLmaxTarget);
			if (ps_time_step_type==CFL_LOCAL) ps_CFLlocal=min(ps_CFLlocal*ps_time_step_ramp_growth,ps_CFLlocalTarget);
		}
	}
	
	double a;
	if (ps_time_step_type==FIXED) {
		for (int gid=0;gid<grid.size();++gid) for (int c=0;c<grid[gid].cellCount;++c) dtau[gid].cell(c)=ps_time_step_current;
	} else if (ps_time_step_type==CFL_MAX) {
		double min_dt=1.e20;
		for (int gid=0;gid<grid.size();++gid) {
			if (equations[gid]==NS) {
				ps_time_step_current=1.e20;
				for (int c=0;c<grid[gid].cellCount;++c) {
					a=ns[gid].material.a(ns[gid].p.cell(c),ns[gid].T.cell(c));
					ps_time_step_current=min(ps_time_step_current,ps_CFLmax*grid[gid].cell[c].lengthScale/(fabs(ns[gid].V.cell(c))+a));
				}
				MPI_Allreduce(&ps_time_step_current,&ps_time_step_current,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
				if (min_dt>ps_time_step_current) min_dt=ps_time_step_current;
			}
			ps_time_step_current=min_dt;
			for (int c=0;c<grid[gid].cellCount;++c) dtau[gid].cell(c)=ps_time_step_current;
		}
	} else if (ps_time_step_type==CFL_LOCAL) {
		// Determine ps_time step with CFL condition
		for (int c=0;c<grid[gid].cellCount;++c) {
			a=ns[gid].material.a(ns[gid].p.cell(c),ns[gid].T.cell(c));
			dtau[gid].cell(c)=ps_CFLlocal*grid[gid].cell[c].lengthScale/(fabs(ns[gid].V.cell(c))+a);
			ps_time_step_current=1.;
		}
	}
	
	// Output current ps_time step, max dt and max CFL
	max_cfl=-1.;
	double min_dt=1.e20;
	if (equations[gid]==NS) {
		for (int c=0;c<grid[gid].cellCount;++c) {
			a=ns[gid].material.a(ns[gid].p.cell(c),ns[gid].T.cell(c));
			max_cfl=max(max_cfl,(fabs(ns[gid].V.cell(c))+a)*ps_time_step_current/grid[gid].cell[c].lengthScale);
			if (ps_time_step_type==CFL_LOCAL) min_dt=min(min_dt,dtau[gid].cell(c));
		}
		if (ps_time_step_type==CFL_LOCAL) {
			MPI_Allreduce(&min_dt,&min_dt,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		} else {
			MPI_Allreduce(&max_cfl,&max_cfl,1, MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		}
	}
	
	if (ps_time_step_type==CFL_LOCAL) {
		ps_time_step_current=min_dt;
		max_cfl=ps_CFLlocal;
	}
	
	return;
}

