/*

 constraints.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2010  Oliver Hahn

 */

 /*
 This file is modified by Jingyuan Chen in 2018.
 */

#include "constraints.hh"

double dsigma2_tophat( double k, void *params );
double dsigma2_gauss( double k, void *params );
double find_coll_z( const std::vector<double>& z, const std::vector<double>& sigma, double nu );
void compute_sigma_tophat( config_file& cf, transfer_function *ptf, double R, std::vector<double>& z, std::vector<double>& sigma );
void compute_sigma_gauss( config_file& cf, transfer_function *ptf, double R, std::vector<double>& z, std::vector<double>& sigma );



double dsigma2_tophat( double k, void *pparams )
{
	if( k<=0.0 )
		return 0.0;

	char **params = reinterpret_cast<char**> (pparams);

	transfer_function *ptf = reinterpret_cast<transfer_function*>(params[0]);
	double x = k * *reinterpret_cast<double*>(params[1]);
	double nspect = *reinterpret_cast<double*>(params[2]);

	double w = 3.0*(sin(x)-x*cos(x))/(x*x*x);

	double tfk = ptf->compute(k,total);
	return k*k * w*w * pow(k,nspect) * tfk*tfk;
}

double dsigma2_gauss( double k, void *pparams )
{
	if( k<=0.0 )
		return 0.0;

	char **params = reinterpret_cast<char**> (pparams);

	transfer_function *ptf = reinterpret_cast<transfer_function*> (params[0]);
	double x = k * *reinterpret_cast<double*>(params[1]);
	double nspect = *reinterpret_cast<double*>(params[2]);

	double w = exp(-x*x*0.5);

	double tfk = ptf->compute(k,total);
	return k*k * w*w * pow(k,nspect) * tfk*tfk;
}

double find_coll_z( const std::vector<double>& z, const std::vector<double>& sigma, double nu )
{
	double dcoll = 1.686/nu;
	double zcoll = 0.0;

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, z.size());

	gsl_spline_init (spline, &sigma[0], &z[0], z.size() );

	zcoll = gsl_spline_eval(spline, dcoll, acc );

	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

	return zcoll;
}


void compute_sigma_tophat( config_file& cf, transfer_function *ptf, double R, std::vector<double>& z, std::vector<double>& sigma )
{
	z.clear();
	sigma.clear();

	cosmology cosm( cf );
	CosmoCalc ccalc( cosm, ptf );

	double zmin = 0.0, zmax = 200.0;
	int nz = 100;
	for( int i=0; i <nz; ++i )
		z.push_back( zmax - i*(zmax-zmin)/(nz-1.0) );

	double D0 = ccalc.CalcGrowthFactor(1.0);

	double sigma8 = cf.getValue<double>("cosmology","sigma_8");
	double nspec = cf.getValue<double>("cosmology","nspec");

	double sigma0 = 0.0;
	{
		double eight=8.0;

		char *params[3];
		params[0] = reinterpret_cast<char*> (ptf);
		params[1] = reinterpret_cast<char*> (&eight);
		params[2] = reinterpret_cast<char*> (&nspec);

		sigma0 = sqrt(4.0*M_PI*integrate( &dsigma2_tophat, 1e-4, 1e4, reinterpret_cast<void*>(params) ));
	}

	for( int i=0; i <nz; ++i )
	{
		void *params[3];
		params[0] = reinterpret_cast<char*> (ptf);
		params[1] = reinterpret_cast<char*> (&R);
		params[2] = reinterpret_cast<char*> (&nspec);

		double sig = sqrt(4.0*M_PI*integrate( &dsigma2_tophat, 1e-4, 1e4, reinterpret_cast<void*>(params) ));
		double Dz  = ccalc.CalcGrowthFactor(1./(1.+z[i]));
		sigma.push_back( sig*sigma8/sigma0*Dz/D0 );
	}
}

void compute_sigma_gauss( config_file& cf, transfer_function *ptf, double R, std::vector<double>& z, std::vector<double>& sigma )
{
	z.clear();
	sigma.clear();

	cosmology cosm( cf );
	CosmoCalc ccalc( cosm, ptf );

	double z_start = cf.getValue<double>("setup","zstart");

	double zmin = z_start, zmax = 200.0;
	int nz = 100;
	for( int i=0; i <nz; ++i )
		z.push_back( zmax - i*(zmax-zmin)/(nz-1.0) );

	double D0 = ccalc.CalcGrowthFactor(1.0);

	double sigma8 = cf.getValue<double>("cosmology","sigma_8");
	double nspec = cf.getValue<double>("cosmology","nspec");

	double sigma0 = 0.0;
	{
		double eight=8.0;

		char *params[3];
		params[0] = reinterpret_cast<char*> (ptf);
		params[1] = reinterpret_cast<char*> (&eight);
		params[2] = reinterpret_cast<char*> (&nspec);

		sigma0 = sqrt(4.0*M_PI*integrate( &dsigma2_tophat, 1e-4, 1e4, reinterpret_cast<void*>(params) ));
	}

	for( int i=0; i <nz; ++i )
	{
		void *params[3];
		params[0] = reinterpret_cast<char*> (ptf);
		params[1] = reinterpret_cast<char*> (&R);
		params[2] = reinterpret_cast<char*> (&nspec);

		double sig = sqrt(4.0*M_PI*integrate( &dsigma2_gauss, 1e-4, 1e4, reinterpret_cast<void*>(params) ));
		double Dz  = ccalc.CalcGrowthFactor(1./(1.+z[i]));

		//std::cerr << z[i] << "    " << sig << std::endl;
		sigma.push_back( sig*sigma8/sigma0*Dz/D0 );
	}
}


constraint_set::constraint_set( config_file& cf, transfer_function *ptf )
: pcf_( &cf ), ptf_( ptf )
{
	pcosmo_ = new Cosmology( cf );
	pccalc_ = new CosmoCalc( *pcosmo_, ptf_ );
	dplus0_ = 1.0;//pccalc_->CalcGrowthFactor( 1.0 );


	unsigned i=0;

    double astart = 1.0/(1.0+pcf_->getValue<double>("setup","zstart"));
    unsigned levelmin = pcf_->getValue<unsigned>("setup","levelmin");
	unsigned levelmin_TF = pcf_->getValueSafe<unsigned>("setup","levelmin_TF",levelmin);
	constr_level_ = pcf_->getValueSafe<unsigned>("constraints","level",levelmin_TF);

	constr_level_ = std::max(constr_level_,levelmin_TF);

	double omegam = pcf_->getValue<double>("cosmology","Omega_m");
	double rhom = omegam*2.77519737e11; //... mean matter density in Msun/Mpc^3

	//... use EdS density for estimation
	//double rhom = 2.77519737e11;


	while(true)
	{
		char temp1[128];
		std::string temp2;
		sprintf(temp1,"constraint[%u].index",i);
		if( cf.containsKey( "constraints", temp1 ) )
		{
			//!... parse a new constraint
			constraint new_c;
			//!... read type index (from 0 to 17)
			new_c.index = cf.getValue<int>( "constraints", temp1 );
			if( new_c.index<0 or new_c.index>17 )
				throw std::runtime_error("Unknown constraint type!\n");


			//!... read position
			sprintf(temp1,"constraint[%u].pos",i);
			temp2 = cf.getValue<std::string>( "constraints", temp1 );
			sscanf(temp2.c_str(), "%lf,%lf,%lf", &new_c.x, &new_c.y, &new_c.z);

			//!... read mass
			sprintf(temp1,"constraint[%u].mass",i);
			double mass = cf.getValue<double>( "constraints", temp1 );
			new_c.Rg = pow((mass/pow(2.*M_PI,1.5)/rhom),1./3.);

			if( new_c.index == 0 )
			{
				//! 0 means peak hight constraint
				//!... peak type constraints take a scale and a peak height

				double Rtophat = pow(mass/4.0*3.0/M_PI/rhom,1./3.);
				sprintf(temp1,"constraint[%u].nu",i);
				double nu = cf.getValue<double>( "constraints", temp1 );

				std::vector<double> z,sigma;
				compute_sigma_tophat( cf, ptf, Rtophat, z, sigma );
				double zcoll = find_coll_z( z, sigma, nu );

				LOGINFO("Probable collapse redshift for constraint %d : z = %f @ M = %g", i, zcoll,mass );

				compute_sigma_gauss( cf, ptf, new_c.Rg, z, sigma );
				new_c.value = nu*sigma.back();

				LOGINFO("Constraint %d : peak with Rg=%g h-1 Mpc and nu = %g",i,new_c.Rg,new_c.value);
				LOGINFO("Constraint %3d : peak",i);
				LOGINFO("   M = %g h-1 M_o, nu = %.2f sigma", mass, nu );
				LOGINFO("   estimated z_coll = %f, sigma = %f", zcoll, new_c.value );

			}
            //! index 1-3 means density gradient constraint. Notice the position should be the same as where the peak is in the following constraints

			if (new_c.index == 1 )
			{
				sprintf(temp1,"constraint[%u].s1",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : s1= %g",i,new_c.value);
			}
			if (new_c.index == 2 )
			{
				sprintf(temp1,"constraint[%u].s2",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : s2= %g",i,new_c.value);
			}
			if (new_c.index == 3 )
			{
				sprintf(temp1,"constraint[%u].s3",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : s3= %g",i,new_c.value);
			}
            //! index 4-9 means Hessian matrix constraint, which control the shape, orientation and compactness of the peak.

			if (new_c.index == 4 )
			{
				sprintf(temp1,"constraint[%u].H11",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : H11= %g",i,new_c.value);
			}
			if (new_c.index == 5 )
			{
				sprintf(temp1,"constraint[%u].H22",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : H22= %g",i,new_c.value);
			}
			if (new_c.index == 6 )
			{
				sprintf(temp1,"constraint[%u].H33",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : H33= %g",i,new_c.value);
			}
			if (new_c.index == 7 )
			{
				sprintf(temp1,"constraint[%u].H12",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : H12= %g",i,new_c.value);
			}
			if (new_c.index == 8 )
			{
				sprintf(temp1,"constraint[%u].H13",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : H13= %g",i,new_c.value);
			}
			if (new_c.index == 9 )
			{
				sprintf(temp1,"constraint[%u].H23",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : H23= %g",i,new_c.value);
			}
            //! index 10-12 means peculiear acceleration constraints.

			if (new_c.index == 10 )
			{
				sprintf(temp1,"constraint[%u].g1",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : g1= %g",i,new_c.value);
			}
			if (new_c.index == 11 )
			{
				sprintf(temp1,"constraint[%u].g2",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : g2= %g",i,new_c.value);
			}
			if (new_c.index == 12 )
			{
				sprintf(temp1,"constraint[%u].g3",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : g3= %g",i,new_c.value);
			}
            //! index 13-17 means tidal matrix constraints.
			if (new_c.index == 13 )
			{
				sprintf(temp1,"constraint[%u].T11",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : T11= %g",i,new_c.value);
			}
			if (new_c.index == 14 )
			{
				sprintf(temp1,"constraint[%u].T22",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : T22= %g",i,new_c.value);
			}
			if (new_c.index == 15 )
			{
				sprintf(temp1,"constraint[%u].T12",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : T12= %g",i,new_c.value);
			}
			if (new_c.index == 16 )
			{
				sprintf(temp1,"constraint[%u].T13",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : T13= %g",i,new_c.value);
			}
			if (new_c.index == 17 )
			{
				sprintf(temp1,"constraint[%u].T23",i);
				new_c.value = cf.getValue<double>( "constraints", temp1 );
				LOGINFO("Constraint %d : T23= %g",i,new_c.value);
			}

			new_c.Rg2 = new_c.Rg*new_c.Rg;

			cset_.push_back( new_c );

		}else
			break;

		++i;
	}

	LOGINFO("Found %d density constraint(s) to be obeyed.",cset_.size());
}


void constraint_set::wnoise_constr_corr( double dx, size_t nx, size_t ny, size_t nz, std::vector<double>& g0, matrix& cinv, fftw_complex* cw )
{
	double lsub = nx*dx;
	double dk = 2.0*M_PI/lsub, d3k=dk*dk*dk;

	double pnorm = pcf_->getValue<double>("cosmology","pnorm");
	double nspec = pcf_->getValue<double>("cosmology","nspec");
	pnorm *= dplus0_*dplus0_;

	size_t nconstr = cset_.size();
	size_t nzp=nz/2+1;



	/*for( size_t i=0; i<nconstr; ++i )
	 for( size_t j=0; j<nconstr; ++j )
	 {
	 std::cerr << "fact    = " << (cset_[j].sigma-g0[j])*cinv(i,j) << "\n";
	 std::cerr << "g(j)    = " << cset_[j].sigma << "\n";
	 std::cerr << "g0(j)   = " << g0[j] << "\n";
	 std::cerr << "qinv    = " << cinv(i,j) << "\n";
	 }
	 */


	double chisq = 0.0, chisq0 = 0.0;
	for( size_t i=0; i<nconstr; ++i )
		for( size_t j=0; j<nconstr; ++j )
		{
			chisq += cset_[i].value*cinv(i,j)*cset_[j].value;
			chisq0 += g0[i]*cinv(i,j)*g0[j];
		}
	LOGINFO("Chi squared for the constraints:\n       sampled = %f, desired = %f", chisq0, chisq );

	std::vector<double> value(nconstr,0.0);

	#pragma omp parallel
	{
		std::vector<double> value_loc(nconstr,0.0);

		#pragma omp for
		for( int ix=0; ix<(int)nx; ++ix )
		{
			double iix(ix); if( iix > nx/2 ) iix-=nx;
			iix *= 2.0*M_PI/nx;

			for( size_t iy=0; iy<ny; ++iy )
			{
				double iiy(iy); if( iiy > ny/2 ) iiy-=ny;
				iiy *= 2.0*M_PI/nx;
				for( size_t iz=0; iz<nzp; ++iz )
				{
					double iiz(iz);
					iiz *= 2.0*M_PI/nx;

					double k = sqrt(iix*iix+iiy*iiy+iiz*iiz)*(double)nx/lsub;

					double T = ptf_->compute(k,total);
					double Pk = pnorm*T*T*pow(k,nspec)*d3k;

					size_t q = ((size_t)ix*ny+(size_t)iy)*nzp+(size_t)iz;

					double fac = sqrt(Pk);

					for( unsigned i=0; i<nconstr; ++i )
						for( unsigned j=0; j<=i; ++j )
						{
							std::complex<double>
							ci = eval_constr(cset_[i].index,i,iix,iiy,iiz),
							cj = eval_constr(cset_[j].index,j,iix,iiy,iiz);

							RE(cw[q]) += (cset_[j].value-g0[j])*cinv(i,j) * std::real(ci)*fac;
							IM(cw[q]) += (cset_[j].value-g0[j])*cinv(i,j) * std::imag(ci)*fac;
							//! Commend the above two lines and use the following lines if you want to generate mean field(under random realization), residual field or true mean field.
							//! See van de Weygaert & Bertschinger 1996 for details.
							//! mean field:
//							RE(cw[q]) = (g0[j])*cinv(i,j) * std::real(ci)*fac;
//							IM(cw[q]) = (g0[j])*cinv(i,j) * std::imag(ci)*fac;

							//! residual field:
//							RE(cw[q]) -= (g0[j])*cinv(i,j) * std::real(ci)*fac;
//							IM(cw[q]) -= (g0[j])*cinv(i,j) * std::imag(ci)*fac;

							//! true mean field:
//							RE(cw[q]) = (cset_[j].value)*cinv(i,j) * std::real(ci)*fac;
//							IM(cw[q]) = (cset_[j].value)*cinv(i,j) * std::imag(ci)*fac;
							if( i!=j )
							{
								RE(cw[q]) += (cset_[i].value-g0[i])*cinv(j,i) * std::real(cj)*fac;
								IM(cw[q]) += (cset_[i].value-g0[i])*cinv(j,i) * std::imag(cj)*fac;
								//! Also commend the above two lines and use the following lines for other field options.
								//! mean field
//								RE(cw[q]) = (g0[i])*cinv(j,i) * std::real(cj)*fac;
//								IM(cw[q]) = (g0[i])*cinv(j,i) * std::imag(cj)*fac;
								//! residual field
//								RE(cw[q]) -= (g0[i])*cinv(j,i) * std::real(cj)*fac;
//								IM(cw[q]) -= (g0[i])*cinv(j,i) * std::imag(cj)*fac;
								//! true mean field
//								RE(cw[q]) = (cset_[i].value)*cinv(j,i) * std::real(cj)*fac;
//								IM(cw[q]) = (cset_[i].value)*cinv(j,i) * std::imag(cj)*fac;

							}
							else
							{
								if( iz>0&&iz<nz/2 )
									value_loc[i] += 2.0*std::real(std::conj(ci)*std::complex<double>(RE(cw[q]),IM(cw[q])))*fac;
								else
									value_loc[i] += std::real(std::conj(ci)*std::complex<double>(RE(cw[q]),IM(cw[q])))*fac;
							}
						}
				}

			}

		}


		//.. 'critical' section for the global reduction
		#pragma omp critical
		{
			for(int i=0; i<(int)nconstr; ++i )
				value[i] += value_loc[i];
		}
	}

	for(int i=0; i<(int)nconstr; ++i )
	{
		switch (cset_[i].index)
		{
			case 0: LOGINFO("Constraint %3d : sigma = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 1: LOGINFO("Constraint %3d : s1 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 2: LOGINFO("Constraint %3d : s2 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 3: LOGINFO("Constraint %3d : s3 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 4: LOGINFO("Constraint %3d : H11 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 5: LOGINFO("Constraint %3d : H22 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 6: LOGINFO("Constraint %3d : H33 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 7: LOGINFO("Constraint %3d : H12 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 8: LOGINFO("Constraint %3d : H13 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 9: LOGINFO("Constraint %3d : H23 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 10: LOGINFO("Constraint %3d : g1 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 11: LOGINFO("Constraint %3d : g2 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 12: LOGINFO("Constraint %3d : g3 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 13: LOGINFO("Constraint %3d : T11 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 14: LOGINFO("Constraint %3d : T22 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 15: LOGINFO("Constraint %3d : T12 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 16: LOGINFO("Constraint %3d : T13 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
			case 17: LOGINFO("Constraint %3d : T23 = %+6f (%+6f)",i,value[i],cset_[i].value); break;
		}
	}
}



void constraint_set::wnoise_constr_corr( double dx, fftw_complex* cw, size_t nx, size_t ny, size_t nz, std::vector<double>& g0 )
{
	size_t nconstr = cset_.size();
	size_t nzp=nz/2+1;

	g0.assign(nconstr,0.0);

	double pnorm = pcf_->getValue<double>("cosmology","pnorm");
	double nspec = pcf_->getValue<double>("cosmology","nspec");
	pnorm *= dplus0_*dplus0_;
	double lsub = nx*dx;
	double dk = 2.0*M_PI/lsub, d3k=dk*dk*dk;
    LOGINFO("Before applying constraints, the values of constraints functions are:");

	for( size_t i=0; i<nconstr; ++i )
	{
		double gg = 0.0;

		#pragma omp parallel for reduction(+:gg)
		for( int ix=0; ix<(int)nx; ++ix )
		{
			double iix(ix); if( iix > nx/2 ) iix-=nx;
			iix *= 2.0*M_PI/nx;

			for( size_t iy=0; iy<ny; ++iy )
			{
				double iiy(iy); if( iiy > ny/2 ) iiy-=ny;
				iiy *= 2.0*M_PI/nx;
				for( size_t iz=0; iz<nzp; ++iz )
				{
					double iiz(iz);
					iiz *= 2.0*M_PI/nx;

					double k = sqrt(iix*iix+iiy*iiy+iiz*iiz)*(double)nx/lsub;
					double T = ptf_->compute(k,total);

					std::complex<double> v(std::conj(eval_constr(cset_[i].index,i,iix,iiy,iiz)));

					v *= sqrt(pnorm*pow(k,nspec)*T*T*d3k);


					if( iz>0&&iz<nz/2)
						v*=2;

					size_t q = ((size_t)ix*ny+(size_t)iy)*nzp+(size_t)iz;

					std::complex<double> ccw(RE(cw[q]),IM(cw[q]));

					gg += std::real(v*ccw);

				}
			}
		}

		g0[i] = gg;
		LOGINFO("g_%d = %.6f",i, g0[i]);
	}
}




void constraint_set::icov_constr( double dx, size_t nx, size_t ny, size_t nz, matrix& cij )
{
	size_t nconstr = cset_.size();
	size_t nzp=nz/2+1;

	double pnorm = pcf_->getValue<double>("cosmology","pnorm");
	double nspec = pcf_->getValue<double>("cosmology","nspec");
	pnorm *= dplus0_*dplus0_;

	cij		= matrix(nconstr,nconstr);

	double lsub = nx*dx;
	double dk = 2.0*M_PI/lsub, d3k=dk*dk*dk;

	//... compute lower triangle of covariance matrix
	//... and fill in upper triangle
	for( unsigned i=0; i<nconstr; ++i )
		for( unsigned j=0; j<=i; ++j )
		{

			float c1(0.0), c2(0.0);

#pragma omp parallel for reduction(+:c1,c2)
			for( int ix=0; ix<(int)nx; ++ix )
			{
				double iix(ix); if( iix > nx/2 ) iix-=nx;
				iix *= 2.0*M_PI/nx;

				for( size_t iy=0; iy<ny; ++iy )
				{
					double iiy(iy); if( iiy > ny/2 ) iiy-=ny;
					iiy *= 2.0*M_PI/nx;
					for( size_t iz=0; iz<nzp; ++iz )
					{
						double iiz(iz);
						iiz *= 2.0*M_PI/nx;

						double k = sqrt(iix*iix+iiy*iiy+iiz*iiz)*(double)nx/lsub;
						double T = ptf_->compute(k,total);
						std::complex<double> v(std::conj(eval_constr(cset_[i].index,i,iix,iiy,iiz)));
						v *= eval_constr(cset_[j].index,j,iix,iiy,iiz);
						v *= pnorm * pow(k,nspec) * T * T * d3k;

						if( iz>0&&iz<nz/2)
							v*=2;

						c1 += std::real(v);
						c2 += std::real(std::conj(v));
					}
				}
			}

			cij(i,j) = c1;
			cij(j,i) = c2;
		}

	//... invert convariance matrix
	cij.invert();

}


void constraint_set::constraint_test( double dx, fftw_complex* cw, size_t nx, size_t ny, size_t nz, std::vector<double>& h0 )
{
	size_t nconstr = cset_.size();
	size_t nzp=nz/2+1;

	h0.assign(18,0.0);

	double pnorm = pcf_->getValue<double>("cosmology","pnorm");
	double nspec = pcf_->getValue<double>("cosmology","nspec");
	pnorm *= dplus0_*dplus0_;
	double lsub = nx*dx;
	double dk = 2.0*M_PI/lsub, d3k=dk*dk*dk;
    LOGINFO("After applying constraints, the values of constraints functions are:");

	for( size_t i=0; i<18; ++i )
	{
		double gg = 0.0;

		#pragma omp parallel for reduction(+:gg)
		for( int ix=0; ix<(int)nx; ++ix )
		{
			double iix(ix); if( iix > nx/2 ) iix-=nx;
			iix *= 2.0*M_PI/nx;

			for( size_t iy=0; iy<ny; ++iy )
			{
				double iiy(iy); if( iiy > ny/2 ) iiy-=ny;
				iiy *= 2.0*M_PI/nx;
				for( size_t iz=0; iz<nzp; ++iz )
				{
					double iiz(iz);
					iiz *= 2.0*M_PI/nx;

					double k = sqrt(iix*iix+iiy*iiy+iiz*iiz)*(double)nx/lsub;
					double T = ptf_->compute(k,total);

					std::complex<double> v(std::conj(eval_constr(i,0,iix,iiy,iiz)));

					v *= sqrt(pnorm*pow(k,nspec)*T*T*d3k);


					if( iz>0&&iz<nz/2)
						v*=2;

					size_t q = ((size_t)ix*ny+(size_t)iy)*nzp+(size_t)iz;

					std::complex<double> ccw(RE(cw[q]),IM(cw[q]));

					gg += std::real(v*ccw);

				}
			}
		}

		h0[i] = gg;
		LOGINFO("g_%d,constr = %.6f",i, h0[i]);
	}
}
