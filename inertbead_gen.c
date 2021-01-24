# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
# include <malloc.h>
# include <time.h>
# include "common.h"
# include <sys/time.h>
# include <mkl_lapacke.h>

#define NMAXCONNECT 5
#define ACCEPTEUR   1
#define DONNEUR     2

#define anint(x) ((x)>0 ? floor((x)+0.5) : ceil((x)-0.5))

double ConvSI(double valeur, char *unite) {

	double pi, boltzmann, avogadro, qElectron, aBohr, Epsilon0, UnSur4PiEpsilon0, debye, planck, planckbar, calorie;

	/*Valeurs des constantes physiques*/
        pi = 4. * atan(1.);
        boltzmann = 1.3806503e-23;
        avogadro = 6.02214199e23;
        qElectron = 1.60217653e-19;
        aBohr = 0.5291772083e-10;
        Epsilon0 = 8.854187e-12;
        UnSur4PiEpsilon0 = 8.98755179e+09;
        debye = 3.33564e-30;
        planck = 6.62606876e-34;
        planckbar = planck / (2. * pi);
        calorie = 4.1868;

	//grandeurs sans unite
	if	(strcmp(unite,"-")		==0)	{return (valeur);}

	//energies -> Joule
	else if	(strcmp(unite,"J")		==0)	{return (valeur);}
	else if	(strcmp(unite,"cal")		==0)	{return (valeur * calorie);}
	else if	(strcmp(unite,"eV")		==0)	{return (valeur * qElectron);}
	else if	(strcmp(unite,"cal/mol")	==0)	{return (valeur * calorie/avogadro);}
	else if	(strcmp(unite,"kcal/mol")	==0)	{return (valeur * 1000*calorie/avogadro);}
	else if	(strcmp(unite,"kJ/mol")		==0)	{return (valeur * 1000/avogadro);}
	else if	(strcmp(unite,"J/mol")		==0)	{return (valeur * 1./avogadro);}
	else if	(strcmp(unite,"K")		==0)	{return (valeur * boltzmann);}

	//Pression -> Pascal
	else if	(strcmp(unite,"Pa")		==0)	{return (valeur);}
	else if	(strcmp(unite,"MPa")		==0)	{return (valeur * 1e6);}
	else if	(strcmp(unite,"GPa")		==0)	{return (valeur * 1e9);}
	else if	(strcmp(unite,"bar")		==0)	{return (valeur * 1e5);}

	//distances -> metre
	else if	(strcmp(unite,"m")		==0)	{return (valeur);}
	else if	(strcmp(unite,"/m")		==0)	{return (valeur);}
	else if	(strcmp(unite,"metre")		==0)	{return (valeur);}
	else if	(strcmp(unite,"/metre")		==0)	{return (valeur);}
	else if	(strcmp(unite,"ang")		==0)	{return (valeur * 1e-10);}
	else if	(strcmp(unite,"/ang")		==0)	{return (valeur * 1e10);}
	else if	(strcmp(unite,"angstrom")	==0)	{return (valeur * 1e-10);}
	else if	(strcmp(unite,"/angstrom")	==0)	{return (valeur * 1e10);}

	//masses -> kg
	else if	(strcmp(unite,"kg")		==0)	{return (valeur);}
	else if	(strcmp(unite,"g")		==0)	{return (valeur * 1e-3);}
	else if	(strcmp(unite,"kg/mol")		==0)	{return (valeur * 1./avogadro);}
	else if	(strcmp(unite,"g/mol")		==0)	{return (valeur * 1e-3/avogadro);}

	//charges electrostatiques -> Coulombs
	else if	(strcmp(unite,"C")		==0)	{return (valeur);}
	else if	(strcmp(unite,"e-")		==0)	{return (valeur * qElectron);}

	//angles -> degres
	else if	(strcmp(unite,"degre")		==0)	{return (valeur);}
	else if	(strcmp(unite,"radian")		==0)	{return (valeur *180./pi);}

	//unites diverses (constantes de force notamment)
	else if	(strcmp(unite,"J/mol.m6")	==0)	{return (valeur * 1./avogadro);}
	else if	(strcmp(unite,"kJ/mol.m6")	==0)	{return (valeur * 1000/avogadro);}
	else if	(strcmp(unite,"J/mol.ang6")	==0)	{return (valeur * 1./avogadro*1e-60);}
	else if	(strcmp(unite,"kJ/mol.ang6")	==0)	{return (valeur * 1000/avogadro*1e-60);}
	else if	(strcmp(unite,"cal/mol.m6")	==0)	{return (valeur * calorie/avogadro);}
	else if	(strcmp(unite,"kcal/mol.m6")	==0)	{return (valeur * 1000*calorie/avogadro);}
	else if	(strcmp(unite,"cal/mol.ang6")	==0)	{return (valeur * calorie/avogadro*1e-60);}
	else if	(strcmp(unite,"kcal/mol.ang6")	==0)	{return (valeur * 1000*calorie/avogadro*1e-60);}
	else if	(strcmp(unite,"kcal/mol/ang")	==0)	{return (valeur * 1000*calorie/avogadro*1e10);}
	else if	(strcmp(unite,"kcal/mol/ang2")	==0)	{return (valeur * 1000*calorie/avogadro*1e20);}
	else if	(strcmp(unite,"kcal/mol/ang3")	==0)	{return (valeur * 1000*calorie/avogadro*1e30);}
	else if	(strcmp(unite,"kcal/mol/ang4")	==0)	{return (valeur * 1000*calorie/avogadro*1e40);}

	//unite non reconnue
	else {
		fprintf(stdout,"Erreur : unite _%s_ non referencee\n",unite);
		exit(1);
	}

}

void readFAtomes(void) {



	FILE	*fatomesFile;
	char	ligne[512], chaine[512], chaine1[512], chaine2[512];

	int   	ligneActive, i, j, numeroType, nombreTermeChampForce;
	int	    nombreTermes, nombreTermeReactivite;
	int	    at1, at2, at3, at4, at5, at6;
	int	    ZmatriceChargee, reactiviteChargee, champDeForcesCharge;
	double	reel, reel1, reel2, reel3, reel4, reel5, reel6;
	char	unite[16], unite1[16], unite2[16], unite3[16], unite4[16];
	
	
	/*Ouverture du fichier*/
	fatomesFile = fopen("FAtomes", "r");
        if (fatomesFile == NULL) {
		fprintf(stderr, "ERREUR : Mauvaise ouverture du fichier %s\n", fatomesFileName);
		exit(1);
        }

	/*Initialisation*/
	ligneActive = 0;
	numeroType  = -1;
	nombreTermes = 0;
	ZmatriceChargee = 0;
	reactiviteChargee = 0;
	champDeForcesCharge = 0;

	/*Lecture du fichier*/
	while (fgets(ligne,512,fatomesFile)) {
	 
		if (strncmp(ligne, "*", 1) == 0) continue;
		if (strncmp(ligne, "#", 1) == 0) continue;

		ligneActive++;

		/*Recuperation du nombre de types d'atomes*/
		if (strncmp(ligne, "NbTypesAtomes", sizeof("NbTypesAtomes") - 1) == 0 && ligneActive == 1) {
			NombreTypesAtomes = atoi(&ligne[sizeof("NbTypesAtomes") - 1]);

			/*Allocation en consequence*/
			nomAtome = (char**)malloc(NombreTypesAtomes*sizeof(char*));
			for ( i = 0; i < NombreTypesAtomes; i++) nomAtome[i] = (char*)malloc(16*sizeof(char));
			masseAtome = (double*)malloc(NombreTypesAtomes*sizeof(double));
			chargeAtome = (double*)malloc(NombreTypesAtomes*sizeof(double));

			longueurLiaison = (double**)malloc(NombreTypesAtomes*sizeof(double*));
			for ( i = 0; i < NombreTypesAtomes; i++) {
				longueurLiaison[i] = (double*)malloc(NombreTypesAtomes*sizeof(double));
				for ( j = 0; j < NombreTypesAtomes; j++) longueurLiaison[i][j] = 0.0;
			}


		} else if (strncmp(ligne, "NbTypesAtomes", sizeof("NbTypesAtomes") - 1) == 0 && ligneActive != 1) {
			fprintf(stdout, "Erreur : Le mot-clef _NbTypesAtomes_ doit etre le premier mot-clef du fichier _FichierAtomes_\n");
			exit(1);
		}

		/*Recuperation du nom du type d'atome*/
		if (strncmp(ligne, "nom", sizeof("nom") - 1) == 0 && strncmp(ligne, "nomXYZ", sizeof("nomXYZ") - 1) != 0) {
			numeroType++;
			sscanf(ligne, "%s %s\n",chaine, nomAtome[numeroType]);
		}
		/*S'il y a un nom XYZ, c'est ce nom que l'on garde*/
		if (strncmp(ligne, "nomXYZ", sizeof("nomXYZ") - 1) == 0) {
			sscanf(ligne, "%s %s\n",chaine, nomAtome[numeroType]);
		}

		/*Verification que le type de la particule est bien Atome*/
		if (strncmp(ligne, "type", sizeof("type") - 1) == 0) {
			sscanf(ligne, "%s %s\n",chaine1, chaine2);
			if (strcmp(chaine2,"Atome") != 0) {
				fprintf(stdout, "Erreur : Toutes les particules doivent etre de type Atome\n");
				exit(1);
			}
		}

		/**Recuperation de la masse du type*/
		if (strncmp(ligne, "masse", sizeof("masse") - 1) == 0) {
			sscanf(ligne, "%s %lf %s\n", chaine, &reel, unite);
			ConvSI(reel,unite);
			masseAtome[numeroType] = reel;
		}

		/*Recuperation de la charge du type*/
		if (strncmp(ligne, "charge", sizeof("charge") - 1) == 0) {
		        sscanf(ligne, "%s %lf %s\n", chaine, &reel, unite);
                        ConvSI(reel,unite);
                        chargeAtome[numeroType] = reel;	
		}

		/*Recuperation de l'ordre des types dans la maille*/
		if (strncmp(ligne, "PositionDesAtomes", sizeof("PositionDesAtomes") - 1) == 0) {
			fgets(ligne, 512, fatomesFile);
			while ((strncmp(ligne, "*", 1) == 0) || (strncmp(ligne, "#", 1) == 0)) fgets(ligne, 512, fatomesFile);
			sscanf(ligne, "%d\n", &nombreAtomeMaille);

			/*Allocation*/
			typeAtomeVsNumeroDsMaille = (int*)malloc(nombreAtomeMaille*sizeof(int));
			reactiviteVsNumeroDsMaille = (int*)malloc(nombreAtomeMaille*sizeof(int));
			for ( i = 0; i< nombreAtomeMaille; i++) reactiviteVsNumeroDsMaille[i] = 0;
            nbAtomesDsMailleDeType = (int*)malloc(NombreTypesAtomes*sizeof(int));
            for ( j = 0; j < NombreTypesAtomes; j++) nbAtomesDsMailleDeType[j] = 0;


			for ( i = 0; i< nombreAtomeMaille; i++) {
				fgets(ligne, 512, fatomesFile);
				sscanf(ligne, "%s %lf %lf %lf\n", chaine, &reel, &reel, &reel);
				for ( j = 0; j < NombreTypesAtomes; j++) {
					if (strncmp(chaine, nomAtome[j], 16) == 0) {
						typeAtomeVsNumeroDsMaille[i] = j;
                        nbAtomesDsMailleDeType[j]++;
						break;
					}
				}
				
			}
		}

		/*Recuperation de la Zmatrice*/
		if (strncmp(ligne, "Zmatrice" , sizeof("Zmatrice")  - 1) == 0) {
			fgets(ligne, 512, fatomesFile);
                        while ((strncmp(ligne, "*", 1) == 0) || (strncmp(ligne, "#", 1) == 0)) fgets(ligne, 512, fatomesFile);
                        sscanf(ligne, "%d\n", &nombreAtomesZmatrice);

			if (nombreAtomesZmatrice != nombreAtomeMaille) {
				fprintf(stdout, "Erreur : le code de posttraitement ne fonctionne que pour nombreAtomesZmatrice (%d) = nombreAtomeMaille (%d) !", nombreAtomesZmatrice, nombreAtomeMaille);
				exit(1);
			}

			/*Allocation*/
			ConnectiviteMaille = (int**) malloc(nombreAtomesZmatrice*sizeof(int*));
			for (i = 0; i < nombreAtomesZmatrice; i++) ConnectiviteMaille[i] = (int*)malloc(NMAXCONNECT*sizeof(int));

			/*Initialisation*/
			for (i = 0; i < nombreAtomesZmatrice; i++) {
				for ( j = 0; j < NMAXCONNECT; j++) ConnectiviteMaille[i][j] = -1;
			}
			

			/*Remplissage de la Zmatrice*/
			for ( i = 0; i < nombreAtomesZmatrice; i++) {
				fgets(ligne, 512, fatomesFile);
				reel1 = reel2 = reel3 = reel4 = reel5 = reel6 = -1;
				sscanf(ligne, "%lf %lf %lf %lf %lf %lf\n",&reel1,&reel2,&reel3,&reel4,&reel5,&reel6);

                                at1 = (int) reel1;
                                at2 = (int) reel2;
                                at3 = (int) reel3;
                                at4 = (int) reel4;
                                at5 = (int) reel5;
                                at6 = (int) reel6;

                                if(at3 == -1) {
                                        ConnectiviteMaille[at1][0] = at2;
                                } else if (at4 == -1) {
                                        ConnectiviteMaille[at1][0] = at2;
                                        ConnectiviteMaille[at1][1] = at3;
                                } else if (at5 == -1){
                                        ConnectiviteMaille[at1][0] = at2;
                                        ConnectiviteMaille[at1][1] = at3;
                                        ConnectiviteMaille[at1][2] = at4;
                                } else if (at6 == -1){
                                        ConnectiviteMaille[at1][0] = at2;
                                        ConnectiviteMaille[at1][1] = at3;
                                        ConnectiviteMaille[at1][2] = at4;
                                        ConnectiviteMaille[at1][3] = at5;
                                } else {
                                        ConnectiviteMaille[at1][0] = at2;
                                        ConnectiviteMaille[at1][1] = at3;
                                        ConnectiviteMaille[at1][2] = at4;
                                        ConnectiviteMaille[at1][3] = at5;
                                        ConnectiviteMaille[at1][4] = at6;
                                }
			}

			ZmatriceChargee = 1;
		}

		/*Recuperation des longueurs de liaison*/
		/* Admettons pour l'instant que le champ de forces est complet*/
		if (strncmp(ligne, "ChampDeForces" , sizeof("ChampDeForces")  - 1) == 0) {

			fgets(ligne, 512, fatomesFile);
                        while ((strncmp(ligne, "*", 1) == 0) || (strncmp(ligne, "#", 1) == 0)) fgets(ligne, 512, fatomesFile);
                        sscanf(ligne, "%d\n", &nombreTermeChampForce);

			nombreTermes = 0;
			while (champDeForcesCharge == 0) {
				fgets(ligne, 512, fatomesFile);
				if (strncmp(ligne, "*", 1) == 0) continue;
				if (strncmp(ligne, "#", 1) == 0) continue;

				/*Mise a jour du nombre de termes lus*/
				nombreTermes++;

				/*Stockage des longueurs de liaisons*/
				if (strncmp(ligne,"harmbond",sizeof("harmbond")-1) == 0) {
					reel1 = 0; reel2 = 0; reel3 = 0; reel4 = 0;
					sscanf(ligne,"%s %s %s %lf %s %lf %s %lf %s %lf %s",
                                                      chaine,chaine1,chaine2,&reel1,unite1,
                                                                             &reel2,unite2,
                                                                             &reel3,unite3,
                                                                             &reel4,unite4);
					
                                	for ( i = 0; i < NombreTypesAtomes; i++) {
                                	        if (strcmp(chaine1, nomAtome[i]) == 0) {
                                	                break;
                                	        }
                                	}
                                	for ( j = 0; j < NombreTypesAtomes; j++) {
                                	        if (strcmp(chaine2, nomAtome[j]) == 0) {
                                	                break;
                                	        }
                                	}
					// conversion en unite SI
					longueurLiaison[i][j] = ConvSI(reel1, unite1);
				//	printf("reel1 and unite1 are:%e %s\n",reel1,unite1);//
					// conversion en Angstrom
					longueurLiaison[i][j] *= 1.0e10;

					longueurLiaison[j][i] = longueurLiaison[i][j];
                //    printf("Bond length is:%d %d %e\n",i,j,longueurLiaison[i][j]);//
				}

				if (nombreTermes == nombreTermeChampForce-1) champDeForcesCharge = 1;
			}
		}

		/*Recuperation de la reactivite*/
		if (strncmp(ligne, "Reactivite" , sizeof("Reactivite")  - 1) == 0) {

                        fgets(ligne, 512, fatomesFile);
                        while ((strncmp(ligne, "*", 1) == 0) || (strncmp(ligne, "#", 1) == 0)) fgets(ligne, 512, fatomesFile);
                        sscanf(ligne, "%d\n", &nombreTermeReactivite);

			if (nombreTermeReactivite != 2) {
				fprintf(stdout, "Erreur : Le nombre de termes dans Reactivite doit etre 2 (pas %d) !\n", nombreTermeReactivite);
				exit(1);
			}

			nombreTermes = 0;
			while (reactiviteChargee == 0) {
				fgets(ligne, 512, fatomesFile);
				if (strncmp(ligne, "*", 1) == 0) continue;
                                if (strncmp(ligne, "#", 1) == 0) continue;

				/*Mise a jour du nombre de termes lus*/
				nombreTermes++;

				reel = -1;
				sscanf(ligne, "%lf %s\n",&reel,chaine);
				at1 = (int) reel;

				if (strncmp(chaine, "ACCEPTEUR", sizeof("ACCEPTEUR")) == 0 ) reactiviteVsNumeroDsMaille[at1] = ACCEPTEUR;
				else if (strncmp(chaine, "DONNEUR", sizeof("DONNEUR")) == 0 ) reactiviteVsNumeroDsMaille[at1] = DONNEUR;

				if (nombreTermes == nombreTermeReactivite) reactiviteChargee = 1;
	
			}

		}


	}
}

void makeArrays(void)
{
int i;//,j,k;

//Atomic data
atometypexyz=(char**)malloc(natomes*sizeof(char*));
for (i=0;i<natomes;i++) atometypexyz[i]=(char*)calloc(3,sizeof(char));

atpos=(double**)malloc(3*sizeof(double*));
for (i=0;i<3;i++) atpos[i]=(double*)calloc(natomes,sizeof(double));

atposreal=(double**)malloc(3*sizeof(double*));
for (i=0;i<3;i++) atposreal[i]=(double*)calloc(natomes,sizeof(double));

//Coordination 
nnb=(int*)calloc(natomes,sizeof(int));
//crd=(double *)calloc(natomes,sizeof(double));

nbtbl=(int **)malloc(NMAXCONNECT*sizeof(int *));
for (i=0;i<NMAXCONNECT;i++) nbtbl[i]=(int *)calloc(natomes,sizeof(int));

//rstor=(double **)malloc(NMAXCONNECT*sizeof(double *));
//for (i=0;i<NMAXCONNECT;i++) rstor[i]=(double *)calloc(natomes,sizeof(double));

//sigstor=(double ***)malloc(3*sizeof(double **));
//for (i=0;i<3;i++) { sigstor[i]=(double **)malloc(NMAXCONNECT*sizeof(double *));
//for (j=0;j<NMAXCONNECT;j++) { sigstor[i][j]=(double *)calloc(natomes,sizeof(double));}}

pbvec=(int**)malloc(3*sizeof(int*));
for (i=0;i<3;i++) pbvec[i]=(int*)calloc(natomes,sizeof(int));

cluster=(int*)calloc(natomes,sizeof(int));
//distcluster=(int*)calloc(natomes,sizeof(int));

}

void freeArrays(void)
{
int i; //,j,k;

//Atomic data
for (i=0;i<3;i++) free(atpos[i]); free(atpos);
for (i=0;i<3;i++) free(atposreal[i]); free(atposreal);
for (i=0;i<3;i++) free(pbvec[i]); free(pbvec);

//molecule data
free(nnb);
for (i=0;i<NMAXCONNECT;i++) free(nbtbl[i]); free(nbtbl);
free(cluster);

}


void complementNghbLists(void) {
	int i,j,k,l;
	
	int ip,jp; //lp;
	int imaille, jmaille;
	int typei, typej;
	
	//FILE *verletxyzFile;
	
	int neighcell[3][14];
	
	int nverlet[3],****verlet;
	int nx,ny,nz;//ntot;
	int icell1,jcell1,kcell1,icell2,jcell2,kcell2;
	int ipcell,jpcell;
	double eps,rverlet0,rverlet[3];
	
	double xi,yi,zi,xj,yj,zj,xij,yij,zij;
	double rsq,r,rcutsq;
	
	//FILE *neighborFile;
	//double xt,yt,zt;
	
	//int icellXmin,icellXmax;
	//int icellYmin,icellYmax;
	//int icellZmin,icellZmax;
	
	//int ix,iy,iz,icell,iat;
	//int jx,jy,jz,jcell,jat;
	//int kx,ky,kz,kcell,kat,kloc;
	//int lx,ly,lz,lcell,lat;
	
	int inbi,inbj;
	int nAtomesVerlet;
	
	//FILE *neighxyzFile;
	
	rcutsq=rcut*rcut;
	
	//#######################################################################
	//##              Allocate each atom to its verlet cell                ##
	//#######################################################################
	
	//minimum radius for verlet cells
	rverlet0 = 0.0;
	for ( i = 0; i < NombreTypesAtomes; i++) {
		for ( j = i; j < NombreTypesAtomes; j++) {
    //printf("Bond length is:%d %d %e \n",i,j,longueurLiaison[i][j]);//
			if (longueurLiaison[i][j] > rverlet0) rverlet0 = longueurLiaison[i][j];
		}
	}
	//rverlet0 += 0.1;
	rverlet0 *= 1.5;
	
	//Calc tot # of Verlet cells
	eps = 1e-6;
	nverlet[0] = (box[0] + eps) / rverlet0;
	nverlet[1] = (box[1] + eps) / rverlet0;
	nverlet[2] = (box[2] + eps) / rverlet0;
	
	rverlet[0] = (box[0] + eps) / nverlet[0];
    //printf("Box size in x and rverlet is:%e %e %e %d\n",box[0],rverlet[0],rverlet0,nverlet[0]);//
	rverlet[1] = (box[1] + eps) / nverlet[1];
	rverlet[2] = (box[2] + eps) / nverlet[2];
	
	
	//printf("\n");
	//printf("\t\ttaille=%g\t%d Verlet cells\tVerlet radius=%g\n",box[0],nverlet[0],rverlet[0]);
	//printf("\t\ttaille=%g\t%d Verlet cells\tVerlet radius=%g\n",box[1],nverlet[1],rverlet[1]);
	//printf("\t\ttaille=%g\t%d Verlet cells\tVerlet radius=%g\n",box[2],nverlet[2],rverlet[2]);
	fprintf(logFile,"\n");
	fprintf(logFile,"\t\ttaille=%g\t%d Verlet cells\tVerlet radius=%g\n",box[0],nverlet[0],rverlet[0]);
	fprintf(logFile,"\t\ttaille=%g\t%d Verlet cells\tVerlet radius=%g\n",box[1],nverlet[1],rverlet[1]);
	fprintf(logFile,"\t\ttaille=%g\t%d Verlet cells\tVerlet radius=%g\n",box[2],nverlet[2],rverlet[2]);
	
	nAtomesVerlet = 20;
	//nAtomesVerlet = 10;
	//Dynamically allocate verlet (4D) - freed at the end of this routine
	verlet = (int****)malloc((nAtomesVerlet+1)*sizeof(int***));	// nAtomesVerlet+1 because index "0" of array stores the tot # of atoms in the Verlet cell : atom index is stored in range [1;nAtomesVerlet]
	for (i = 0; i < (nAtomesVerlet+1); i++)	{verlet[i]      = (int***)malloc(nverlet[0]*sizeof(int**)); 
	for (j = 0; j < nverlet[0]; j++)     	{verlet[i][j]   = (int**) malloc(nverlet[1]*sizeof(int*)); 
	for (k = 0; k < nverlet[1]; k++)	{verlet[i][j][k]= (int*)  calloc(nverlet[2],sizeof(int));
	}}}
	
	//Initialise to empty Verlet cells
	for (i=0;i<nverlet[0];i++){for (j=0;j<nverlet[1];j++){for (k=0;k<nverlet[2];k++){verlet[0][i][j][k]=0;}}}	
	
	for (ip=0; ip<natomes; ip++){
       // printf("Position of atom %d and verlet r is:%lf %lf\n",ip,atpos[0][ip],rverlet[0]);//	
		nx = atpos[0][ip]/rverlet[0];	//Calc coordinates of Verlet cell for atom ip
        //printf("nx is:%d\n",nx);//
		ny = atpos[1][ip]/rverlet[1];
		nz = atpos[2][ip]/rverlet[2];
   //     printf("Atom number is:%d\n",ip);//	
		//Test on the coordinates of the Verlet cell
		if (nx<0 || nx>nverlet[0]) {
			printf("ERROR in VERLET for atom #%d : nx=%d not in allowed range [%d;%d] calculated from system boundaries for this proc\n",ip,nx,0,nverlet[0]-1);
			fprintf(logFile,"ERROR in VERLET for atom #%d : nx=%d not in allowed range [%d;%d] calculated from system boundaries for this proc\n",ip,nx,0,nverlet[0]-1);
			return;
		}
		if (ny<0 || ny>nverlet[1]) {
			printf("ERROR in VERLET for atom #%d : nx=%d not in allowed range [%d;%d] calculated from system boundaries for this proc\n",ip,ny,0,nverlet[1]-1);
			fprintf(logFile,"ERROR in VERLET for atom #%d : nx=%d not in allowed range [%d;%d] calculated from system boundaries for this proc\n",ip,ny,0,nverlet[1]-1);
			return;
		}
		if (nz<0 || nz>nverlet[2]) {
			printf("ERROR in VERLET for atom #%d : nx=%d not in allowed range [%d;%d] calculated from system boundaries for this proc\n",ip,nz,0,nverlet[2]-1);
			fprintf(logFile,"ERROR in VERLET for atom #%d : nx=%d not in allowed range [%d;%d] calculated from system boundaries for this proc\n",ip,nz,0,nverlet[2]-1);
			return;
		}
	
		verlet[0][nx][ny][nz]++;				//Increment # of atoms in Verlet cell
		if (verlet[0][nx][ny][nz] > nAtomesVerlet) {
			printf("ERROR in makeLists : # of atoms allocated to Verlet cell [%d;%d;%d] = %d",nx,ny,nz,verlet[0][nx][ny][nz]);
			printf(" is larger than nAtomesVerlet = %d \n",nAtomesVerlet);
			printf("\t-> increase nAtomesVerlet in main.h\n");
			fprintf(logFile,"ERROR in makeLists : # of atoms allocated to Verlet cell [%d;%d;%d] = %d",nx,ny,nz,verlet[0][nx][ny][nz]);
			fprintf(logFile," is larger than nAtomesVerlet = %d \n",nAtomesVerlet);
			fprintf(logFile,"\t-> increase nAtomesVerlet in main.h\n");
			return;
		}
		verlet[verlet[0][nx][ny][nz]][nx][ny][nz]=ip;	//Store atom number in Verlet cell
	}
	
	#ifdef INFO
	//Print occupancy of each Verlet cell
		fprintf(logFile,"%4d * %4d * %4d verlet cells\n",nverlet[0],nverlet[1],nverlet[2]);
		ntot=0;
		for (i=0;i<nverlet[0];i++){
		for (j=0;j<nverlet[1];j++){
		for (k=0;k<nverlet[2];k++){
		fprintf(logFile,"VerletCell[%4d;%4d;%4d] contains %4d atoms\n",i,j,k,verlet[0][i][j][k]);
		}}}
	#endif
	
	#ifdef INFO
	//Print all Verlet cells separately in .xyz file
		verletxyzFile=fopen("verlet.xyz","w");
		ntot=0;
		for (i=0;i<nverlet[0];i++){
		for (j=0;j<nverlet[1];j++){
		for (k=0;k<nverlet[2];k++){
			ntot+=verlet[0][i][j][k];
			fprintf(verletxyzFile,"%d\t\n%d %d %d\n",verlet[0][i][j][k],i,j,k);
			for (l=1;l<=verlet[0][i][j][k];l++){
				lp=verlet[l][i][j][k];
				fprintf(verletxyzFile,"C\t%lf %lf %lf\t%d\n",atpos[0][lp],atpos[1][lp],atpos[2][lp],verlet[l][i][j][k]);
			}
		}}}
		printf("Ntot=%d\n",ntot);
		fclose(verletxyzFile);
	#endif
	
	//#######################################################################
	//##       For each atom pair, calculate interatomic distance r.       ##
	//#######################################################################
	
	//List of relative coordinates of half the neighbouring cells of a single Verlet cell (including Verlet cell itself)
	//----------------------current cell---------------------------
	neighcell[0][0]  = 0;   neighcell[1][0] =  0;    neighcell[2][0] = 0;
	//------------------------1st layer----------------------------
	neighcell[0][1]  = 0;  neighcell[1][1]  =  0;   neighcell[2][1]  =  1;
	neighcell[0][2]  = 0;  neighcell[1][2]  =  1;   neighcell[2][2]  = -1;
	neighcell[0][3]  = 0;  neighcell[1][3]  =  1;   neighcell[2][3]  =  0;
	neighcell[0][4]  = 0;  neighcell[1][4]  =  1;   neighcell[2][4]  =  1;
	neighcell[0][5]  = 1;  neighcell[1][5]  = -1;   neighcell[2][5]  = -1;
	neighcell[0][6]  = 1;  neighcell[1][6]  = -1;   neighcell[2][6]  =  0;
	neighcell[0][7]  = 1;  neighcell[1][7]  = -1;   neighcell[2][7]  =  1;
	neighcell[0][8]  = 1;  neighcell[1][8]  =  0;   neighcell[2][8]  = -1;
	neighcell[0][9]  = 1;  neighcell[1][9]  =  0;   neighcell[2][9]  =  0;
	neighcell[0][10] = 1;  neighcell[1][10] =  0;   neighcell[2][10] =  1;
	neighcell[0][11] = 1;  neighcell[1][11] =  1;   neighcell[2][11] = -1;
	neighcell[0][12] = 1;  neighcell[1][12] =  1;   neighcell[2][12] =  0;
	neighcell[0][13] = 1;  neighcell[1][13] =  1;   neighcell[2][13] =  1;
	
	//Initialise some tables
	//Initialise the neighbor table nbtbl from ConnectiviteMaille
	k = 0;
   // printf("nMaille is:%d\n",nMailles);
	for (i = 0; i < nMailles; i++ ) {
		for (j = 0; j < nombreAtomeMaille; j++) {
			nnb[k] = -1;
			for (l = 0; l < NMAXCONNECT; l++) {
	///nbtblprintf("Connectivite for %d and %d is : %d\n",j,l,ConnectiviteMaille[j][l]); //
				if (ConnectiviteMaille[j][l] == -1) nbtbl[l][k] = ConnectiviteMaille[j][l];
				else nbtbl[l][k] = i*nombreAtomeMaille + ConnectiviteMaille[j][l];
			}
//	printf("Table of neighbours of  %d and %d is : %d\n",l,k,nbtbl[l][k]); //
			for (l = 0; l < NMAXCONNECT; l++) {
				if (ConnectiviteMaille[j][l] == -1) {
					nnb[k] = l;
					break;
				}
			}
			if (nnb[k] == -1) nnb[k] = NMAXCONNECT;
//	printf("Number of neighbors of atom %d is:%d\n",k,nnb[k]); //
			k++;
		}
	}
	//for (i = 0; i < natomes; i++) crd[i] = 0;
	
	//Start running over Verlet cells
	for (icell1 = 0; icell1 < nverlet[0]; icell1++){	//First (triple-)loop over all Verlet cells on the current proc
	for (jcell1 = 0; jcell1 < nverlet[1]; jcell1++){
	for (kcell1 = 0; kcell1 < nverlet[2]; kcell1++){
	
		for (l=0;l<14;l++) {			//Second loop over half the neighbours (14-63-172 for 1-2-3 layers respectively) of the current Verlet cell - only half to avoid double-counting of pair interactions
	
			icell2=icell1+neighcell[0][l];	//Neighbouring cell coordinates
			jcell2=jcell1+neighcell[1][l];
			kcell2=kcell1+neighcell[2][l];
	
			//Apply PBC to Verlet cells
			pbc[0]=pbc[1]=pbc[2]=0.;
	
			if (npbx==1) {
				if (icell2<0) 			{icell2+=nverlet[0]; pbc[0]=-box[0]; }
				else if (icell2>=nverlet[0])	{icell2-=nverlet[0]; pbc[0]=box[0]; }
			} else {
				if (icell2<0 || icell2>=nverlet[0]) continue;
			}
	
			if (npby==1) {
				if (jcell2<0) 			{jcell2+=nverlet[1]; pbc[1]=-box[1]; }
				else if (jcell2>=nverlet[1])	{jcell2-=nverlet[1]; pbc[1]=box[1]; }
			} else {
				if (jcell2<0 || jcell2>=nverlet[1]) continue;
			}
	
			if (npbz==1) {
				if (kcell2<0)			{kcell2+=nverlet[2]; pbc[2]=-box[2]; }
				else if (kcell2>=nverlet[2])	{kcell2-=nverlet[2]; pbc[2]=box[2]; }
			} else {
				if (kcell2<0 || kcell2>=nverlet[2]) continue;
			}
	
			for (ipcell=1;ipcell<=verlet[0][icell1][jcell1][kcell1];ipcell++){	//Loop over all the atoms in the first Verlet cell
	
				ip=verlet[ipcell][icell1][jcell1][kcell1];
	
				xi=atpos[0][ip];
				yi=atpos[1][ip];
				zi=atpos[2][ip];
	
				for (jpcell=1;jpcell<=verlet[0][icell2][jcell2][kcell2];jpcell++){	//Loop over all the atoms in the neighbouring cell
	
					jp=verlet[jpcell][icell2][jcell2][kcell2];
	
					if (l==0 && ip>=jp) continue;					//Avoids self-interaction and double-counting of pair interactions inside the first Verlet cell

					imaille = ip % nombreAtomeMaille;
					jmaille = jp % nombreAtomeMaille;
					if ((reactiviteVsNumeroDsMaille[imaille] == ACCEPTEUR && reactiviteVsNumeroDsMaille[jmaille] == DONNEUR) || 
					    (reactiviteVsNumeroDsMaille[imaille] == DONNEUR   && reactiviteVsNumeroDsMaille[jmaille] == ACCEPTEUR) ) {
	
						xj = atpos[0][jp] + pbc[0];	//Apply PBC
						yj = atpos[1][jp] + pbc[1];
						zj = atpos[2][jp] + pbc[2];
	
						rsq = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi);		// calculate r_{ij}^2
                        printf("Square root is:%lf\n",rsq);	
						typei = typeAtomeVsNumeroDsMaille[imaille];
						typej = typeAtomeVsNumeroDsMaille[jmaille];
						rcut  = longueurLiaison[typei][typej]*1.5; 
						rcutsq=rcut*rcut;
	
						if (rsq<rcutsq) {
	
							r = sqrt(rsq);
	
							xij=(xj-xi)/r;					//Unit vector from ip to jp
							yij=(yj-yi)/r;					//
							zij=(zj-zi)/r;					//
	
							nnb[ip]++;					//	for ip
							inbi = nnb[ip]-1;				//
							nbtbl[inbi][ip] = jp;				//	fill neigh list if :	-set_fulllist_construct is set (filling both central and ghost atom lists)
							//rstor[inbi][ip]=r;				//				-set_fulllist_construct is not set and ip is a central atom
							//sigstor[0][inbi][ip]=xij;			//
							//sigstor[1][inbi][ip]=yij;			//
							//sigstor[2][inbi][ip]=zij;			//
	
							nnb[jp]++;					//	for jp
                      //      printf("Number of neighbors of atom %d is:%d\n",ip,nnb[ip]); //
							inbj = nnb[jp]-1;				//
							nbtbl[inbj][jp] = ip;				//	fill neigh list if :	-set_fulllist_construct is set (filling both central and ghost atom lists)
							//rstor[inbj][jp]=r;				//				-set_fulllist_construct is not set and jp is a central atom
							//sigstor[0][inbj][jp]=-xij;			//
							//sigstor[1][inbj][jp]=-yij;			//
							//sigstor[2][inbj][jp]=-zij;			//
						}
					}
				}
			}
		}
	}}}
	
	
	//Free verlet (4D)
	for (i=0;i<(nAtomesVerlet+1);i++) {
		for (j=0;j<nverlet[0];j++) {
			for (k=0;k<nverlet[1];k++) {
				free(verlet[i][j][k]);
			}
			free(verlet[i][j]);
		}
		free(verlet[i]);
	}
	free(verlet);
	
	
	//Miscellaneous checks and print outs ##
	
	#ifdef INFO
	//Print coordination dependent XYZ file of all processor atoms
		neighxyzFile=fopen("neighXYZ.xyz","a");
		fprintf(neighxyzFile,"%d\n%lf %lf %lf\n",natomes,box[0],box[1],box[2]);
		for (iat=0;iat<natomesCell;iat++) {
			if (nnb[iat]==0) fprintf(neighxyzFile,"Ne\t%lf %lf %lf\n",atpos[0][iat],atpos[1][iat],atpos[2][iat]);
			if (nnb[iat]==1) fprintf(neighxyzFile,"F \t%lf %lf %lf\n",atpos[0][iat],atpos[1][iat],atpos[2][iat]);
			if (nnb[iat]==2) fprintf(neighxyzFile,"O \t%lf %lf %lf\n",atpos[0][iat],atpos[1][iat],atpos[2][iat]);
			if (nnb[iat]==3) fprintf(neighxyzFile,"N \t%lf %lf %lf\n",atpos[0][iat],atpos[1][iat],atpos[2][iat]);
			if (nnb[iat]==4) fprintf(neighxyzFile,"C \t%lf %lf %lf\n",atpos[0][iat],atpos[1][iat],atpos[2][iat]);
			if (nnb[iat]> 4) fprintf(neighxyzFile,"S \t%lf %lf %lf\n",atpos[0][iat],atpos[1][iat],atpos[2][iat]);
		}
		fclose(neighxyzFile);
	#endif
	
	#ifdef INFO
	//Print neigh lists.
		printf("\n");
		printf("#########################\n");
		printf("##Print neighbour lists##\n");
		printf("#########################\n");
		printf("\n");
	
		for (i=0;i<natomes;i++) {
			printf("%5d (%s)\t%d voisins :\n",i,atometypexyz[i],nnb[i]);
			for (j=1;j<=nnb[i];j++) printf("\t%8d%8d(%s)\t%e\n",j,nbtbl[j][i],atometypexyz[nbtbl[j][i]],rstor[j][i]);
		}
		printf("\n");
	#endif
	
	#ifdef INFO
	//Write XYZ file containing SR-MR neighbours around atom "i"
		neighborFile=fopen("neighbour.xyz","w");
		i=0;
		fprintf(neighborFile,"%d\n",1+nnb[i]);
		fprintf(neighborFile,"neighbours\n");
		fprintf(neighborFile,"O\t%lf\t%lf\t%lf\n",atpos[0][i],atpos[1][i],atpos[2][i]);
		for (j=0;j<nnb[i];j++) {
			jat=nbtbl[j][i];
			xt=atpos[0][jat];
			yt=atpos[1][jat];
			zt=atpos[2][jat];
			fprintf(neighborFile,"C\t%lf\t%lf\t%lf\n",xt,yt,zt);
		}
		fclose(neighborFile);
	#endif
	
	
}

void searchNeighbors(int iat, int icluster)
{
int j,jat;
//printf("Number of neighbors of atom %d is:%d\n",iat,nnb[iat]); 
for (j = 0; j < nnb[iat]; j++) {
	jat = nbtbl[j][iat];
   // printf("Neighbour's table is:%d\n",nbtbl[j][iat]);
	if (cluster[jat] >= 0) continue;
	cluster[jat] = icluster;

	pbvec[0][jat] = pbvec[0][iat] - (int)((atpos[0][jat]-atpos[0][iat])/(0.5*box[0]));
	pbvec[1][jat] = pbvec[1][iat] - (int)((atpos[1][jat]-atpos[1][iat])/(0.5*box[1]));
	pbvec[2][jat] = pbvec[2][iat] - (int)((atpos[2][jat]-atpos[2][iat])/(0.5*box[2]));

	searchNeighbors(jat,icluster);
}
//printf("Pbvectors for atom %d are:%lf %lf %lf\n",iat,pbvec[0][iat],pbvec[1][iat],pbvec[2][iat]); 
}


main()
{
	int i,j,k,l,nbeadstot;
	int nbeadspmol,lamda;
	int natpm;
	int ii,jj;
	int *molid_bead,ibead;
    	int iatomes;
	int imol,nmols;
	int iat,maxsteps,z;
    	int ictheta,ictheta_i,ictheta_j,nthetamax;
    	int idr,nrmax;
    	int  shiftX,shiftY,shiftZ;
	char **id;
	char fName[34],outName[34];
    	char **type;
	double rc,rc2,dr;
    	double ctheta,ctheta_i,ctheta_j,dctheta;
	double *mass_bead,*invmass_bead;
    	double ****kount_avg,****kount_intra_avg,****kount_12_avg,****kount_13_avg;
	double mass1,mass2,mass3;
	double inertia[3][3];
    	double **fpos;
	double rx_ib[2400],ry_ib[2400],rz_ib[2400];
	double fx_ib[2400],fy_ib[2400],fz_ib[2400];
	double rx_CG[12][200],ry_CG[12][200],rz_CG[12][200];
    	double fx_CG[12][200],fy_CG[12][200],fz_CG[12][200];
	double drx,dry,drz;
    	double principalDirx[12][200],principalDiry[12][200],principalDirz[12][200];
    	double principalx[2400],principaly[2400],principalz[2400];
	double *work;
	double matrice[3][3];
	double eigenvalue[3];
	double inertia1,inertia2,inertia3;
	double exijpar,eyijpar,ezijpar;
	double exijper,eyijper;
	double pspar,psper;
	double *massat;
	double ****fmpar_avg,****fmper_avg,****fmpar_12_avg,****fmper_12_avg,****fmpar_13_avg,****fmper_13_avg,****fmpar_intra_avg,****fmper_intra_avg;
	double rij,rij2;
	double rxij,ryij,rzij;
	FILE *infile,*outfile,*uncoorfile;
	lapack_int lwork; 	

    	z=0;
    	maxsteps=0;
    /* Reading  FAtomes file */
    	npbx=npby=npbz=1;
    	readFAtomes();
	rc=20.0;
	rc2=rc*rc;
	dr=0.2;
	dctheta=0.2;
    	nrmax=(int)(rc/dr);
    	nthetamax=(int)(2.0/dctheta);
		fmpar_avg=(double ****)malloc(nrmax*sizeof(double ***));
		fmper_avg=(double ****)malloc(nrmax*sizeof(double ***));
		fmpar_12_avg=(double ****)malloc(nrmax*sizeof(double ***));
		fmpar_13_avg=(double ****)malloc(nrmax*sizeof(double ***));
		fmpar_intra_avg=(double ****)malloc(nrmax*sizeof(double ***));
		fmper_12_avg=(double ****)malloc(nrmax*sizeof(double ***));
		fmper_13_avg=(double ****)malloc(nrmax*sizeof(double ***));
		fmper_intra_avg=(double ****)malloc(nrmax*sizeof(double ***));
		kount_avg=(double ****)malloc(nrmax*sizeof(double ***));
		kount_intra_avg=(double ****)malloc(nrmax*sizeof(double ***));
		kount_12_avg=(double ****)malloc(nrmax*sizeof(double ***));
		kount_13_avg=(double ****)malloc(nrmax*sizeof(double ***));
		for (i=0;i<nrmax;i++) {
		fmpar_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		fmpar_12_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		fmpar_13_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		fmpar_intra_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		fmper_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		fmper_12_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		fmper_13_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		fmper_intra_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		kount_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		kount_intra_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		kount_12_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		kount_13_avg[i]=(double ***)malloc(nthetamax*sizeof(double **));
		for (j=0;j<nthetamax;j++) {
		fmpar_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		fmpar_12_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		fmpar_13_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		fmpar_intra_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		fmper_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		fmper_12_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		fmper_13_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		fmper_intra_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		kount_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		kount_intra_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		kount_12_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		kount_13_avg[i][j]=(double **)malloc(nthetamax*sizeof(double *));
		for (k=0;k<nthetamax;k++) {
		fmpar_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		fmpar_12_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		fmpar_13_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		fmpar_intra_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		fmper_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		fmper_12_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		fmper_13_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		fmper_intra_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		kount_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		kount_intra_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		kount_12_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
		kount_13_avg[i][j][k]=(double *)malloc(nthetamax*sizeof(double ));
					}  
				}  
			}  
	printf("Give maximum number of steps:\n");
	scanf("%d",&maxsteps);
//	printf("Give no. of molecules:\n");
//	scanf("%d",&nmols);
//	printf("Give total number of beads per molecules:\n");
//	scanf("%d",&nbeadspmol);
/*	printf("Spatial resolution dr (in Angstroms) = %f\n", dr);
	printf("Spatial Range rc (in Angstroms)  = %f\n", rc); */
	printf("Give coarse-graining level:\n");
	scanf("%d",&lamda);
//	uncoorfile=fopen("uncoor.dat","w");
	
	for (z=0;z<maxsteps;z=z+1) {
		sprintf(fName,"PasDeCalcul_Iteration_%d.xyz",z);
		infile=fopen(fName,"r");	
		if (infile==NULL)
		{
			printf("Can't open %s\n",fName);
			exit(1);
		} 
	
	printf("Reading file %s\n",fName);	
	//printf("Total number of atoms are:%d\n",natomes);
	uncoorfile=fopen("uncoor.dat","w");
    	lwork=3*3-1;
    	work=(double*)malloc(lwork*sizeof(double)); 
	inertia1=0.0;
	inertia2=0.0;
	inertia3=0.0;
    
//	mass1=15.034e-03;
//	mass2=13.019e-03;
//	mass3=14.027e-03;
    	logFile=fopen("printMolecules.log","w");

		iat=0;
            atposMIN[0] = atposMIN[1] = atposMIN[2] = 1e8;
            atposMAX[0] = atposMAX[1] = atposMAX[2] = -1e8;
	 //printf("Inertia tensor element  at the start of frame %d is:%lf\n",z,inertia[0][0]);//
	fscanf(infile,"%d\n",&natomes);
	//printf("Total number of atoms are:%d\n",natomes);
    	nMailles=(int) (natomes/nombreAtomeMaille);
	fscanf(infile,"%lf %lf %lf\n",&box[0],&box[1],&box[2]);
    	makeArrays();
	massat=(double *)malloc(natomes*sizeof(double ));
	type=(char **)malloc(natomes*sizeof(char *));
	for (i=0;i<natomes;i++)
	type[i]=(char *)malloc(7*sizeof(char));
    	id=(char **)malloc(natomes*sizeof(char *));
    	for (i=0;i<natomes;i++) id[i]=(char *)malloc(7*sizeof(char)); 
    	fpos=(double **)malloc(3*sizeof(double *));
    	for (i=0;i<3;i++)
    	fpos[i]=(double *)malloc(natomes*sizeof(double ));
	j=-1;
	/*Read atomic positions */
    	for (iatomes=0;iatomes < natomes;iatomes++) {
	fscanf(infile,"%s %lf %lf %lf %s",type[iatomes],&atposreal[0][iat],&atposreal[1][iat],&atposreal[2][iat],id[iatomes]);
	for (j=0;j<NombreTypesAtomes;j++) {
	if (strcmp(type[iatomes],nomAtome[j])==0) massat[iatomes]=masseAtome[j];
	}
//	printf("Atom id is:%d, atom mass is:%lf\n",iatomes,massat[iatomes]);
    	atpos[0][iat]=atposreal[0][iat];
    	atpos[1][iat]=atposreal[1][iat];
    	atpos[2][iat]=atposreal[2][iat];
      	if (atpos[0][iat] < atposMIN[0]) {atposMIN[0] = atpos[0][iat];
                } else if (atpos[0][iat] > atposMAX[0]) { atposMAX[0] = atpos[0][iat];								}

                if (atpos[1][iat] < atposMIN[1]) {atposMIN[1] = atpos[1][iat];
                } else if (atpos[1][iat] > atposMAX[1]) { atposMAX[1] = atpos[1][iat];								}

                if (atpos[2][iat] < atposMIN[2]) {atposMIN[2] = atpos[2][iat];
                } else if (atpos[2][iat] > atposMAX[2]) { atposMAX[2] = atpos[2][iat];								}

    
    		iat++;
		}
            //Manage PBC or nonPBC (bound system in range [0;atposMAX-atposMIN] along each direction)
            for (iat = 0; iat < natomes; iat++) {
            	//PBC along X
            	if (npbx == 1) {
            		if (atposreal[0][iat] < 0) {shiftX = atposreal[0][iat]/box[0]-1;}
            		else if (atposreal[0][iat] >= box[0]) {shiftX = atposreal[0][iat]/box[0];}
            		else {shiftX = 0;}
            
            		atpos[0][iat] -= shiftX*box[0];
            
            		if (atpos[0][iat] < 0 || atpos[0][iat] >= box[0]) printf("PBC not working!!!\n");
            	} else {
            		atpos[0][iat] -= atposMIN[0];
            	}
            
            	//PBC along Y
            	if (npby==1) {
            		if (atposreal[1][iat] < 0) {shiftY = atposreal[1][iat]/box[1]-1;}
            		else if (atposreal[1][iat] >= box[1]) {shiftY = atposreal[1][iat]/box[1];}
            		else {shiftY = 0;}
            
            		atpos[1][iat] -= shiftY*box[1];
            
            		if (atpos[1][iat] < 0 || atpos[1][iat] >= box[1]) printf("PBC not working!!!\n");
            	} else {
            		atpos[1][iat] -= atposMIN[1];
           	}
            
            	//PBC along Z
            	if (npbz==1) {
            		if (atposreal[2][iat] < 0) {shiftZ = atposreal[2][iat]/box[2]-1;}
            		else if (atposreal[2][iat] >= box[2]) {shiftZ = atposreal[2][iat]/box[2];}
            		else {shiftZ = 0;}
            
            		atpos[2][iat] -= shiftZ*box[2];
            
            		if (atpos[2][iat] < 0 || atpos[2][iat] >= box[2]) printf("PBC not working!!!\n");
            	} else {
            		atpos[2][iat] -= atposMIN[2];
            	}
            }
    //Build neighbor lists//
    fprintf(logFile,"\n\tbuilding neighbour lists\n");
    complementNghbLists();


    //Construction of the molecules//
    icluster=-1;
    for (iat=0;iat<natomes;iat++) {cluster[iat]=-1;}
    
    for (iat=0;iat<natomes;iat++) {
        if (cluster[iat]>=0) continue;
        pbvec[0][iat]=pbvec[1][iat]=pbvec[2][iat]=0;
        icluster++;
        cluster[iat]=icluster;
        searchNeighbors(iat,icluster);
        }
        icluster++;
	nmols=icluster;
	printf("Total number of molecules is:%d\n",icluster);
	printf("Total number of atom types is:%d\n",NombreTypesAtomes);
	for (j=0;j<NombreTypesAtomes;j++) 
	printf("Atomic masses are:%lf\n",masseAtome[j]);
	natpm=natomes/nmols;
	nbeadspmol=natpm/lamda;	
	nbeadstot=nbeadspmol*nmols;

    	molid_bead=(int *)malloc(nbeadstot*sizeof(int)); 
    	mass_bead=(double *)malloc(nbeadspmol*sizeof(double));
    	invmass_bead=(double *)malloc(nbeadspmol*sizeof(double));

	for (ii=0;ii<nbeadspmol;ii++) mass_bead[ii]=0.0;
	for (ii=0;ii<nbeadspmol;ii++) { 
	    for (imol=0;imol<nmols;imol++) {
            fx_CG[ii][imol]=0.0;
            fy_CG[ii][imol]=0.0;
            fz_CG[ii][imol]=0.0;
	        rx_CG[ii][imol]=0.0; 
	        ry_CG[ii][imol]=0.0; 
	        rz_CG[ii][imol]=0.0; 
            principalDirx[ii][imol]=0.0;
            principalDiry[ii][imol]=0.0;
            principalDirz[ii][imol]=0.0;
		}
	}

	for (i=0;i<nbeadstot;i++) {
		fx_ib[i]=0.0;
		fy_ib[i]=0.0;
		fz_ib[i]=0.0;
		rx_ib[i]=0.0;
		ry_ib[i]=0.0;
		rz_ib[i]=0.0;
		principalx[i]=0.0;
		principaly[i]=0.0;
		principalz[i]=0.0;
		}
    	imol=-1;
    	ibead=-1;
     	for (i=0;i<nmols;i++) {
          imol=imol+1;
           for (j=0;j<nbeadspmol;j++) {
                ibead=ibead+1;
                molid_bead[ibead]=imol;
                }
        }	
/*	jj=0;
	for (i=0;i<nbeadspmol;i++) {
		for (j=0;j<lamda;j++) {
		jj=jj+1;
		mass_bead[i]=mass_bead[i]+massat  */     
printf("Status 0 ok.\n");
    for (iat=0;iat<natomes;iat++) {
        atpos[0][iat]+=npbx*pbvec[0][iat]*box[0]; 
        atpos[1][iat]+=npby*pbvec[1][iat]*box[1]; 
        atpos[2][iat]+=npbz*pbvec[2][iat]*box[2]; 
        fprintf(uncoorfile,"%6s %11.8lf %11.8lf %11.8lf\n",type[iat],atpos[0][iat],atpos[1][iat],atpos[2][iat]);
        }
	jj=-1;
	printf("Status 1 ok.\n");
	/* Loop over all molecules */
	for (imol=0;imol<nmols;imol++) {
	// Loop over all beads constituting the chain //
		for (ii=0;ii<nbeadspmol;ii++) {
		invmass_bead[ii]=0.0;
		mass_bead[ii]=0.0;
        // Do not take into account the previous beads //
        fx_CG[ii][imol]=0.0;
        fy_CG[ii][imol]=0.0;
        fz_CG[ii][imol]=0.0;
		rx_CG[ii][imol]=0.0;
		ry_CG[ii][imol]=0.0;
		rz_CG[ii][imol]=0.0;
		 //Loup over the atoms of each bead //
		for (iat=0;iat<lamda;iat++) {
		jj=jj+1;
//		if (strcmp(type[jj],"CH3sp3")==0) {
		mass_bead[ii]=mass_bead[ii]+massat[jj];
		rx_CG[ii][imol]=rx_CG[ii][imol]+massat[jj]*atpos[0][jj];
		ry_CG[ii][imol]=ry_CG[ii][imol]+massat[jj]*atpos[1][jj];
		rz_CG[ii][imol]=rz_CG[ii][imol]+massat[jj]*atpos[2][jj];
	//	}
/*		else if (strcmp(type[jj],"CHsp2")==0) {
		mass_bead[ii]=mass_bead[ii]+mass2;
		rx_CG[ii][imol]=rx_CG[ii][imol]+mass2*atpos[0][jj];
		ry_CG[ii][imol]=ry_CG[ii][imol]+mass2*atpos[1][jj];
		rz_CG[ii][imol]=rz_CG[ii][imol]+mass2*atpos[2][jj];
		}
		else  {
		mass_bead[ii]=mass_bead[ii]+mass3;
		rx_CG[ii][imol]=rx_CG[ii][imol]+mass3*atpos[0][jj];
		ry_CG[ii][imol]=ry_CG[ii][imol]+mass3*atpos[1][jj];
		rz_CG[ii][imol]=rz_CG[ii][imol]+mass3*atpos[2][jj];
		}  
    		fx_CG[ii][imol]=fx_CG[ii][imol]+fpos[0][jj];
    		fy_CG[ii][imol]=fy_CG[ii][imol]+fpos[1][jj];
    		fz_CG[ii][imol]=fz_CG[ii][imol]+fpos[2][jj];  */
    
		}
		invmass_bead[ii]=1.00/mass_bead[ii];
		rx_CG[ii][imol]=rx_CG[ii][imol]*invmass_bead[ii];
		ry_CG[ii][imol]=ry_CG[ii][imol]*invmass_bead[ii];
		rz_CG[ii][imol]=rz_CG[ii][imol]*invmass_bead[ii];
	//	printf("Mass of bead %d is:%lf\n",ii,mass_bead[ii]);
		}
	}
 	/* Calculate inertia tensor for each bead*/ 
	ii=0;
	iat=0;
	jj=-1;
	fclose(uncoorfile);
	for (imol=0;imol<nmols;imol++) {
		for (ii=0;ii<nbeadspmol;ii++) {
			for (i=0;i<3;i++){
				for (j=0;j<3;j++){
			// Zero values for inertia elements before each new bead //
				inertia[i][j]=0.0;
					}
						}
		// Loop over the atoms constituting each bead //
			for (iat=0;iat<lamda;iat++) {
			jj=jj+1;
			drx=atpos[0][jj]-rx_CG[ii][imol];
			dry=atpos[1][jj]-ry_CG[ii][imol];
			drz=atpos[2][jj]-rz_CG[ii][imol];
		 // Fill in the inertia tensor //
//	if (strcmp(type[jj],"CH3sp3")==0) {
	inertia[0][0]=inertia[0][0]+massat[jj]*dry*dry+massat[jj]*drz*drz;
	inertia[1][1]=inertia[1][1]+massat[jj]*drx*drx+massat[jj]*drz*drz;
	inertia[2][2]=inertia[2][2]+massat[jj]*drx*drx+massat[jj]*dry*dry;
	inertia[1][2]=inertia[1][2]-massat[jj]*dry*drz;
	inertia[0][2]=inertia[0][2]-massat[jj]*drx*drz;
	inertia[0][1]=inertia[0][1]-massat[jj]*drx*dry;
/*	}
	else if (strcmp(type[jj],"CHsp2")==0) {
	inertia[0][0]=inertia[0][0]+mass2*dry*dry+mass2*drz*drz;
	inertia[1][1]=inertia[1][1]+mass2*drx*drx+mass2*drz*drz;
	inertia[2][2]=inertia[2][2]+mass2*drx*drx+mass2*dry*dry;
	inertia[1][2]=inertia[1][2]-mass2*dry*drz;
	inertia[0][1]=inertia[0][1]-mass2*drx*dry;
	inertia[0][2]=inertia[0][2]-mass2*drx*drz;
	}
	else {
	inertia[0][0]=inertia[0][0]+mass3*dry*dry+mass3*drz*drz;
	inertia[1][1]=inertia[1][1]+mass3*drx*drx+mass3*drz*drz;
	inertia[2][2]=inertia[2][2]+mass3*drx*drx+mass3*dry*dry;
	inertia[1][2]=inertia[1][2]-mass3*dry*drz;
	inertia[0][2]=inertia[0][2]-mass3*drx*drz;
	inertia[0][1]=inertia[0][1]-mass3*drx*dry;
	} */
}
	inertia[2][1]=inertia[1][2];
	inertia[2][0]=inertia[0][2];
	inertia[1][0]=inertia[0][1];
	 // Diagonalise inertia tensor of the current bead //
	for (i=0;i<3;i++){
		for (j=0;j<3;j++) {
		matrice[i][j]=inertia[i][j];
			}
		}  
	LAPACKE_dsyev_work(LAPACK_ROW_MAJOR, 'V', 'U', 3, *matrice,3, eigenvalue, work, lwork);


    	principalDirx[ii][imol]=matrice[0][0];
    	principalDiry[ii][imol]=matrice[1][0];
    	principalDirz[ii][imol]=matrice[2][0];
		}
	} 
	fclose(infile);
	freeArrays();
	free(work);
	printf("Eigenvalues of the inertia tensor are:%lf %lf %lf\n",eigenvalue[0],eigenvalue[1],eigenvalue[2]);
	for (i=0;i<natomes;i++) {
		free(id[i]);
		free(type[i]);
			}
	free(id);
	free(type);
	for (i=0;i<3;i++) free(fpos[i]);
	free(fpos);
	free(mass_bead);
	free(invmass_bead);	
    	ibead=-1;	

        for (imol=0;imol<nmols;imol++) {
    		for (ii=0;ii<nbeadspmol;ii++) {
                ibead=ibead+1;
                rx_ib[ibead]=rx_CG[ii][imol];
                ry_ib[ibead]=ry_CG[ii][imol]; 
                rz_ib[ibead]=rz_CG[ii][imol];
                fx_ib[ibead]=fx_CG[ii][imol];
                fy_ib[ibead]=fy_CG[ii][imol]; 
                fz_ib[ibead]=fz_CG[ii][imol];
                principalx[ibead]=principalDirx[ii][imol];
                principaly[ibead]=principalDiry[ii][imol];
                principalz[ibead]=principalDirz[ii][imol];
                }
	    }
    	for (i=0;i<nbeadstot-1;i++) {
       		 for (j=i+1;j<nbeadstot;j++){
            // Calculate the distance rij between the CG beads //
            rxij=rx_ib[i]-rx_ib[j];
            ryij=ry_ib[i]-ry_ib[j];
            rzij=rz_ib[i]-rz_ib[j];
            rxij=rxij-box[0]*anint(rxij/box[0]);
            ryij=ryij-box[1]*anint(ryij/box[1]);
            rzij=rzij-box[2]*anint(rzij/box[2]);
            rij2=rxij*rxij+ryij*ryij+rzij*rzij;
            if (rij2 < rc2) {
            rij=sqrt(rij2);
            exijpar=rxij/rij;
            eyijpar=ryij/rij;
            ezijpar=rzij/rij;
            // Projection of the force to eij // 
            exijper=-eyijpar/(sqrt(exijpar*exijpar+eyijpar*eyijpar));
            eyijper= exijpar/(sqrt(exijpar*exijpar+eyijpar*eyijpar));
            
            pspar=(fx_ib[i]-fx_ib[j])*exijpar+(fy_ib[i]-fy_ib[j])*eyijpar+(fz_ib[i]-fz_ib[j])*ezijpar;
            psper=(fx_ib[i]-fx_ib[j])*exijper+(fy_ib[i]-fy_ib[j])*eyijper;
            idr= (int)(rij/dr);
            ctheta=principalx[i]*principalx[j]+principaly[i]*principaly[j]+principalz[i]*principalz[j];

            ctheta_i=-(principalx[i]*exijpar+principaly[i]*eyijpar+principalz[i]*ezijpar);
            ctheta_j=principalx[j]*exijpar+principaly[j]*eyijpar+principalz[j]*ezijpar;
            ictheta=(int) ((ctheta+1.0)/dctheta); 
            ictheta_i=(int) ((ctheta_i+1.0)/dctheta); 
            ictheta_j=(int) ((ctheta_j+1.0)/dctheta); 
		
            if (molid_bead[i]==molid_bead[j])  {
                if (j==(i+1)) {
                // Intramolecular 1-2 interaction //
                    fmpar_12_avg[idr][ictheta_i][ictheta_j][ictheta]=fmpar_12_avg[idr][ictheta_i][ictheta_j][ictheta]+pspar;
                    fmper_12_avg[idr][ictheta_i][ictheta_j][ictheta]=fmper_12_avg[idr][ictheta_i][ictheta_j][ictheta]+psper; 
                    kount_12_avg[idr][ictheta_i][ictheta_j][ictheta]=kount_12_avg[idr][ictheta_i][ictheta_j][ictheta]+1.0;
                    }
                else if (j==(i+2)) {
                    fmpar_13_avg[idr][ictheta_i][ictheta_j][ictheta]=fmpar_13_avg[idr][ictheta_i][ictheta_j][ictheta]+pspar;
                    fmper_13_avg[idr][ictheta_i][ictheta_j][ictheta]=fmper_13_avg[idr][ictheta_i][ictheta_j][ictheta]+psper; 
                    kount_13_avg[idr][ictheta_i][ictheta_j][ictheta]=kount_13_avg[idr][ictheta_i][ictheta_j][ictheta]+1.0;
                    }
                else {
                    fmpar_intra_avg[idr][ictheta_i][ictheta_j][ictheta]=fmpar_intra_avg[idr][ictheta_i][ictheta_j][ictheta]+pspar;
                    fmper_intra_avg[idr][ictheta_i][ictheta_j][ictheta]=fmper_intra_avg[idr][ictheta_i][ictheta_j][ictheta]+psper; 
                    kount_intra_avg[idr][ictheta_i][ictheta_j][ictheta]=kount_intra_avg[idr][ictheta_i][ictheta_j][ictheta]+1.0;
                    }
		}
            else {
                    fmpar_avg[idr][ictheta_i][ictheta_j][ictheta]=fmpar_avg[idr][ictheta_i][ictheta_j][ictheta]+pspar;
                    fmper_avg[idr][ictheta_i][ictheta_j][ictheta]=fmper_avg[idr][ictheta_i][ictheta_j][ictheta]+psper; 
                    kount_avg[idr][ictheta_i][ictheta_j][ictheta]=kount_avg[idr][ictheta_i][ictheta_j][ictheta]+1.0;
                    }
               }  
            }
        }
		
/*End of the analysis for the current frame */ 
}
printf("Status 2 ok.\n");
	/* Store data from all frames */
            sprintf(outName,"MCF_3th_%d.out",z);
            outfile=fopen(outName,"w");
	if (outfile==NULL) 
		{
			printf("Can't open outfile \n");
			exit(1);
		} 
            printf("Writing conservative force components in file MCF_3theta.out\n");
            fprintf(outfile,"Dr=%.2lfAng, Nrmax=%d, dctheta=%lf, Nthetamax=%d\n", dr, nrmax, dctheta, nthetamax);
            fprintf(outfile," \n");
            fprintf(outfile," Average force from inter- and intramolecular contributions\n");
            fprintf(outfile," column: Finter_par Finter_per kount_iner Fintra_par Fintra_per kount_intra\n");
            fprintf(outfile," for a given R value, all costheta, then next R...\n");
            fprintf(outfile," Newton \n");

            for (i=0;i<nrmax;i++) {
                for (j=0;j<nthetamax;j++) {
               		 for (k=0;k<nthetamax;k++) {
               		 	for (l=0;l<nthetamax;l++) {
                    fprintf(outfile,"\t%4d \t%4d \t%4d \t%4d \t%e \t%e \t%10.1lf \t%e \t%e \t%10.1lf \t%e  \t%e \t%10.1lf \t%e \t%e \t%10.1lf\n",i,j,k,l,0.5*fmpar_avg[i][j][k][l],0.5*fmper_avg[i][j][k][l],kount_avg[i][j][k][l],0.5*fmpar_intra_avg[i][j][k][l], 0.5*fmper_intra_avg[i][j][k][l],kount_intra_avg[i][j][k][l],0.5*fmpar_12_avg[i][j][k][l],0.5*fmper_12_avg[i][j][k][l],kount_12_avg[i][j][k][l],0.5*fmpar_13_avg[i][j][k][l],0.5*fmper_13_avg[i][j][k][l],kount_13_avg[i][j][k][l]);  
                                  			   }
                               		           }  
                    			  }  
                         	  }  
fclose(outfile);
	for (i=0;i<nrmax;i++) {
		for (j=0;j<nthetamax;j++) {
			for (k=0;k<nthetamax;k++) {
free(fmpar_avg[i][j][k]);
free(fmper_avg[i][j][k]);
free(fmpar_intra_avg[i][j][k]);
free(fmpar_12_avg[i][j][k]);
free(fmpar_13_avg[i][j][k]);
free(fmper_intra_avg[i][j][k]);
free(fmper_12_avg[i][j][k]);
free(fmper_13_avg[i][j][k]);
free(kount_avg[i][j][k]);
free(kount_intra_avg[i][j][k]);
free(kount_12_avg[i][j][k]);
free(kount_13_avg[i][j][k]);
}
free(fmpar_avg[i][j]);
free(fmpar_intra_avg[i][j]);
free(fmpar_12_avg[i][j]);
free(fmpar_13_avg[i][j]);
free(fmper_avg[i][j]);
free(fmper_intra_avg[i][j]);
free(fmper_12_avg[i][j]);
free(fmper_13_avg[i][j]);
free(kount_avg[i][j]);
free(kount_intra_avg[i][j]);
free(kount_12_avg[i][j]);
free(kount_13_avg[i][j]);
}
free(fmpar_avg[i]);
free(fmpar_intra_avg[i]);
free(fmpar_12_avg[i]);
free(fmpar_13_avg[i]);
free(fmper_avg[i]);
free(fmper_intra_avg[i]);
free(fmper_12_avg[i]);
free(fmper_13_avg[i]);
free(kount_avg[i]);
free(kount_intra_avg[i]);
free(kount_12_avg[i]);
free(kount_13_avg[i]);
}
free(fmpar_avg);
free(fmpar_intra_avg);
free(fmpar_12_avg);
free(fmpar_13_avg);
free(fmper_avg);
free(fmper_intra_avg);
free(fmper_12_avg);
free(fmper_13_avg);
free(kount_avg);
free(kount_intra_avg);
free(kount_12_avg);
free(kount_13_avg);
}
