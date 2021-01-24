# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
# include <malloc.h>
# include <time.h>
# include "common.h"
# include <sys/time.h>
# include <omp.h>

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
for (i=0;i<natomes;i++) free(atometypexyz[i]); free(atometypexyz);
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


int main(int argc, char *argv[])
{
	int i,j,k,nbeadstot;
	int nbeadspmol,lamda;
	int ii,jj;
	int l;
	int *molid_bead,ibead;
    int iatomes;
	int imol,nmols;
	int iat,maxsteps,z;
    int ictheta,nbinctheta;
    int nbonds,nbends;
    int idr,nbin;
    int iconf,cbin,bin;
    int  shiftX,shiftY,shiftZ;
	char idat[4],**id;
	char fName[34],grbondName[34],grbendName[34],grintraName[34],grinterName[34],grtotName[34];
    char ligne[512];
    char **type;
	char atype[8];
	double rmax,dr,vShell;
	double n1,n2;
	double costheta;
    double ctheta,dctheta;
	double *mass_bead,*invmass_bead;
	double mass1,mass2,mass3;
    double **fpos;
	double rx_ib[2400],ry_ib[2400],rz_ib[2400];
	double rx_CG[12][200],ry_CG[12][200],rz_CG[12][200];
	double **gDeRtot_thread,**gDeRinter_thread,**gDeRintra_thread,*gDeRbond,*gDeRbend,*gDeRinter,*gDeRintra,*gDeRtot;
    double norme;
	double dist;
	double xrij,yrij,zrij;
	double ri,rij,rij2;
	double rxij,ryij,rzij,rxij2,ryij2,rzij2;
	FILE *infile,*grbondfile,*grbendfile,*grintrafile,*grinterfile,*grtotfile,*beadfile;
    int tid,nthreads,ji,i8,j8;
   long long int ij,npairs;
    
    /* Recuperation du nombre de threads */
    nthreads=1;
    #pragma omp parallel
    {
    #pragma omp master
    {
    nthreads=omp_get_num_threads();
    }
    }
    	z=0;
    	maxsteps=0;
    /* Reading of FAtomes file */
	npbx=npby=npbz=1;
  /*liste des arguments d'entree */
	while (++i<argc) {
	printf("\n");
        if (strcmp(argv[i],"-dmaxs")==0)
        maxsteps=atoi(argv[++i]);
        else if (strcmp(argv[i],"-dmol")==0)
        nmols=atoi(argv[++i]);
        else if (strcmp(argv[i],"-dbeads")==0)
        nbeadspmol=atoi(argv[++i]);
        else if (strcmp(argv[i],"-dlamda")==0)
        lamda=atoi(argv[++i]);
        else
        printf("Unknown options\n");
        }

    /* Reading of FAtomes file */
    readFAtomes();
	dr=0.02;
	dctheta=0.01;
/*	printf("Give maximum number of steps:\n");
	scanf("%d",&maxsteps);
//	maxsteps=10000000;
	printf("Give no. of molecules:\n");
	scanf("%d",&nmols);
//	nmols=200;
	printf("Give total number of beads per molecules:\n");
	scanf("%d",&nbeadspmol);
//	nbeadspmol=12;
	printf("Give coarse-graining level:\n");
	scanf("%d",&lamda);
//	lamda=4; */
    	iconf=0;

        z=100;
		sprintf(fName,"PasDeCalcul_Iteration_%d.xyz",z);
		infile=fopen(fName,"r");	
		if (infile==NULL)
		{
			printf("Can't open %s\n",fName);
			exit(1);
		} 
	
  fgets(ligne,512,infile);
  sscanf(ligne,"%lli\n",&natomes);
  fgets(ligne,512,infile);
  sscanf(ligne,"%lf%lf%lf\n",&box[0],&box[1],&box[2]);
  fclose(infile);

    rmax=box[0]/2.0;
//  nbpts=(int) (rmax/dr);
    nbin=(int) (rmax/dr)+1;
    nbinctheta=(int) (2.0/dctheta)+1;
//    npbx=npby=npbz=1;
    // Allocation of gDer
  /* gDeRinter */
    gDeRinter = (double*)malloc(nbin*sizeof(double));
    gDeRinter_thread = (double**)malloc(nbin*sizeof(double*));
            for (bin = 0; bin < nbin; bin++){
             gDeRinter_thread[bin] = (double*)malloc(nthreads*sizeof(double));;
            for (j = 0; j < nthreads; j++)  {
            gDeRinter_thread[bin][j] = 0.0;;
                }
        }
  /* gDeRtot */
  /*  gDeRtot = (double*)malloc(nbin*sizeof(double));
    gDeRtot_thread = (double**)malloc(nbin*sizeof(double*));
            for (bin = 0; bin < nbin; bin++){
             gDeRtot_thread[bin] = (double*)malloc(nthreads*sizeof(double));;
            for (j = 0; j < nthreads; j++){
             gDeRtot_thread[bin][j] = 0.0;;
                    }
            }  */
    /* gDeRintra */
    gDeRintra = (double*)malloc(nbin*sizeof(double));
    gDeRintra_thread = (double**)malloc(nbin*sizeof(double*));
            for (bin = 0; bin < nbin; bin++) {
             gDeRintra_thread[bin] = (double*)malloc(nthreads*sizeof(double));;
            for (j = 0; j < nthreads; j++) {
             gDeRintra_thread[bin][j] = 0.0;;
                }
            }
   /* gDeRbond */
    gDeRbond = (double*)malloc(nbin*sizeof(double));

   /* gDeCthetaBend */
    gDeRbend = (double*)malloc(nbinctheta*sizeof(double));

    z=100; 
	for (z=100;z<maxsteps;z=z+100) {
		sprintf(fName,"PasDeCalcul_Iteration_%d.xyz",z);
        	infile=fopen(fName,"r"); 
		    printf("Reading file %s\n",fName);	
  
    nbeadstot=nmols*nbeadspmol;
    molid_bead=(int *)malloc(nbeadstot*sizeof(int)); 
    mass_bead=(double *)malloc(nbeadspmol*sizeof(double));
    invmass_bead=(double *)malloc(nbeadspmol*sizeof(double));
	mass1=15.034e-03;
	mass2=13.019e-03;
	mass3=14.027e-03;
    logFile=fopen("printMolecules.log","w");
	for (ii=0;ii<nbeadspmol;ii++) mass_bead[ii]=0.0;
	for (ii=0;ii<nbeadspmol;ii++) { 
	    for (imol=0;imol<nmols;imol++) {
	        rx_CG[ii][imol]=0.0; 
	        ry_CG[ii][imol]=0.0; 
	        rz_CG[ii][imol]=0.0; 
		}
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
	for (i=0;i<nbeadstot;i++) {
		rx_ib[i]=0.0;
		ry_ib[i]=0.0;
		rz_ib[i]=0.0;
		}
	iat=0;
            atposMIN[0] = atposMIN[1] = atposMIN[2] = 1e8;
            atposMAX[0] = atposMAX[1] = atposMAX[2] = -1e8;
	 //printf("Inertia tensor element  at the start of frame %d is:%lf\n",z,inertia[0][0]);//
	fscanf(infile,"%d\n",&natomes);
	//printf("Total number of atoms are:%d\n",natomes);
    nMailles=(int) (natomes/nombreAtomeMaille);
	fscanf(infile,"%lf %lf %lf\n",&box[0],&box[1],&box[2]);
    makeArrays();
//    nbpts=(int) (rmax/dr);

	type=(char **)malloc(natomes*sizeof(char *));
	for (i=0;i<natomes;i++)
	type[i]=(char *)malloc(7*sizeof(char));
    id=(char **)malloc(natomes*sizeof(char *));
    for (i=0;i<natomes;i++) id[i]=(char *)malloc(7*sizeof(char)); 
    fpos=(double **)malloc(3*sizeof(double *));
    for (i=0;i<3;i++)
    fpos[i]=(double *)malloc(natomes*sizeof(double ));
	//printf("Status 1 ok\n");
	/*Read atomic positions */
    for (iatomes=0;iatomes < natomes;iatomes++) {
	fscanf(infile,"%s %lf %lf %lf %lf %lf %lf %s",type[iatomes],&atposreal[0][iat],&atposreal[1][iat],&atposreal[2][iat],&fpos[0][iat],&fpos[1][iat],&fpos[2][iat],id[iatomes]);
    atpos[0][iat]=atposreal[0][iat];
    atpos[1][iat]=atposreal[1][iat];
    atpos[2][iat]=atposreal[2][iat];
      if (atpos[0][iat] < atposMIN[0]) {atposMIN[0] = atpos[0][iat];
                } else if (atpos[0][iat] > atposMAX[0]) { atposMAX[0] = atpos[0][iat];}

                if (atpos[1][iat] < atposMIN[1]) {atposMIN[1] = atpos[1][iat];
                } else if (atpos[1][iat] > atposMAX[1]) { atposMAX[1] = atpos[1][iat];}

                if (atpos[2][iat] < atposMIN[2]) {atposMIN[2] = atpos[2][iat];
                } else if (atpos[2][iat] > atposMAX[2]) { atposMAX[2] = atpos[2][iat];}

    
   // printf("Atomic position %4d is:%lf %lf %lf\n",iat,atpos[0][iat],atpos[1][iat],atpos[2][iat]);
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

    for (iat=0;iat<natomes;iat++) {
        atpos[0][iat]+=npbx*pbvec[0][iat]*box[0]; 
        atpos[1][iat]+=npby*pbvec[1][iat]*box[1]; 
        atpos[2][iat]+=npbz*pbvec[2][iat]*box[2]; 
        }  
	printf("Status 1 ok\n");
    jj=-1;
	/* Loop over all molecules */
	for (imol=0;imol<nmols;imol++) {
	// Loop over all beads constituting the chain //
		for (ii=0;ii<nbeadspmol;ii++) {
		invmass_bead[ii]=0.0;
		mass_bead[ii]=0.0;
        // Do not take into account the previous beads //
		rx_CG[ii][imol]=0.0;
		ry_CG[ii][imol]=0.0;
		rz_CG[ii][imol]=0.0;
		 //Loup over the atoms of each bead //
		for (iat=0;iat<lamda;iat++) {
        jj=jj+1;
//	strcpy(type[iat],atype);
	if (strcmp(type[jj],"CH3sp3")==0) {
	mass_bead[ii]=mass_bead[ii]+mass1;
	rx_CG[ii][imol]=rx_CG[ii][imol]+mass1*atpos[0][jj];
	ry_CG[ii][imol]=ry_CG[ii][imol]+mass1*atpos[1][jj];
	rz_CG[ii][imol]=rz_CG[ii][imol]+mass1*atpos[2][jj];
	}
	else if (strcmp(type[jj],"CHsp2")==0) {
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
}
	invmass_bead[ii]=1.00/mass_bead[ii];
	rx_CG[ii][imol]=rx_CG[ii][imol]*invmass_bead[ii];
	ry_CG[ii][imol]=ry_CG[ii][imol]*invmass_bead[ii];
	rz_CG[ii][imol]=rz_CG[ii][imol]*invmass_bead[ii];
	}
}
 
//	fclose(uncoorfile);
	fclose(infile);
	freeArrays();
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
//	beadfile=fopen("bead.dat","w");
        for (imol=0;imol<nmols;imol++) {
    		for (ii=0;ii<nbeadspmol;ii++) {
                ibead=ibead+1;
                rx_ib[ibead]=rx_CG[ii][imol];
                ry_ib[ibead]=ry_CG[ii][imol]; 
                rz_ib[ibead]=rz_CG[ii][imol];
//		fprintf(beadfile,"%lf %lf %lf\n",rx_ib[ibead],ry_ib[ibead],rz_ib[ibead]);
                }
    }
 //  fclose(beadfile); 
    // Calculate grbond
    ii=0;
    for (i=0;i<nmols;i++) {
        for (j=0;j<(nbeadspmol-1);j++){

            rxij=rx_ib[ii]-rx_ib[ii+1];
            ryij=ry_ib[ii]-ry_ib[ii+1];
            rzij=rz_ib[ii]-rz_ib[ii+1];
            rij2=rxij*rxij+ryij*ryij+rzij*rzij;
            rij=sqrt(rij2);
//	printf("Bonded beads' distance is:%lf\n",rij);
            // contributions bonding
            if (rij < (rmax+0.0000001)) {
            nbonds=nbonds+1;
            idr=(int)(rij/dr);
//	printf("Status 4  ok\n");
            gDeRbond[idr]=gDeRbond[idr]+1.0;
            }
            ii=ii+1;
            }
            ii=ii+1;
            }
            // calculate grbend // 
	
//	printf("Status 3  ok\n");
            ii=0; 
            for (i=0;i<nmols;i++) {
                  for (j=0;j<(nbeadspmol-2);j++) {
                     rxij=rx_ib[ii]-rx_ib[ii+1];
                     ryij=ry_ib[ii]-ry_ib[ii+1];
                     rzij=rz_ib[ii]-rz_ib[ii+1];
                     rxij2=rx_ib[ii+2]-rx_ib[ii+1];
                     ryij2=ry_ib[ii+2]-ry_ib[ii+1];
                     rzij2=rz_ib[ii+2]-rz_ib[ii+1];
		//	printf("Triplet of beads is %d %d %d\n",ii,ii+1,ii+2);
                     n1=sqrt((rxij*rxij)+(ryij*ryij)+(rzij*rzij));
                     n2=sqrt((rxij2*rxij2)+(ryij2*ryij2)+(rzij2*rzij2));
    
                     costheta=((rxij*rxij2)+(ryij*ryij2)+(rzij*rzij2))/(n1*n2);
                    //  printf("Costheta by bead triplet is:%lf\n",costheta);
                      // contributions bending
                       nbends=nbends+1;
                      ictheta=(int) ((costheta+1.0)/dctheta); 
                      gDeRbend[ictheta]=gDeRbend[ictheta]+1.0;
                      ii=ii+1;
                                              }
                      ii=ii+2;
                                   }  // imol
	
//		printf("Status  ok\n");
            /* Loop ends for current frame */ 
            // Calculate nonbonding g(r)-Double loop over atoms
            npairs=(nbeadstot*(nbeadstot-1))/2;
//		printf("Npairs is:%d\n",npairs);
            #pragma omp parallel private(tid,ij,ji,k,l,bin,xrij,yrij,zrij,dist,i8,j8)
            {

    nthreads=omp_get_num_threads();
	printf("Running on %d threads\n",nthreads);
            /*Recuperation du numero de la thread*/
            tid=omp_get_thread_num();

            #pragma omp for 
            for (ij=0;ij<npairs;ij++) {
		
            /*variables privees */
          /*  long long int ji,i8,j8;
	    double xrij,yrij,zrij;
            double xi,yi,zi;
	    int ix,iy,iz;
            int k,l,bin;
            double dist; */
	
	
	//printf("Status 3  ok\n");
            /*Recuperation des numeros atomiques 'a la JMT' */
            ji=npairs-1-ij;
//		printf("Total no. of pairs is:%d\n",npairs);
            if (ji==0) {
                i8=1;
                j8=0;
            } else {
                i8=floor(sqrt(2*ji-0.25)+0.5);

                if ((i8*(i8+1))/2<=ji) {
                i8=i8+1;
                } else if (ji<(i8*(i8-1))/2) {
                i8=i8-1;
                }
                j8=ji-(i8*(i8-1))/2;
                }
                k=(int) (nbeadstot-i8-1);
                l=(int) (nbeadstot-j8-1);

                xrij=rx_ib[k]-rx_ib[l];
                yrij=ry_ib[k]-ry_ib[l];
                zrij=rz_ib[k]-rz_ib[l];

        	xrij=xrij-box[0]*anint(xrij/box[0]);
		yrij=yrij-box[1]*anint(yrij/box[1]);
		zrij=zrij-box[2]*anint(zrij/box[2]);  

				dist=((xrij*xrij)+(yrij*yrij)+(zrij*zrij));
				dist=sqrt(dist);

      /*          xrij=xrij-box[0]*anint(xrij/box[0]);
                  yrij=yrij-box[1]*anint(yrij/box[1]);
                  zrij=zrij-box[2]*anint(zrij/box[2]);   */ 

//	printf("Current pair of beads is %d and %d\n",k,l);
           //     distance=boxmax;
           /*    for (ix=-npbx;ix<npbx;ix++) {
                    for (iy=-npby;iy<=npby;iy++) {
                        for (iz=-npbz;iz<=npbz;iz++) {
                            xj=rx_ib[jat]+ix*box[0];
                            yj=ry_ib[jat]+iy*box[1];
                            xj=rz_ib[jat]+iz*box[2];
                            tmp=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj);
                            tmp=sqrt(tmp);
                            if (tmp<distance) distance=tmp;
                        }
                    }
                }   */
                // intra+intermolecular contributions

                // intra-molecular contribution
                if (molid_bead[k]==molid_bead[l]) { 
                   if (abs(k-l) > 2) {
                        if (dist<(rmax+0.0000001)) {
                       bin=(int) (dist/dr);
                        gDeRintra_thread[bin][tid]=gDeRintra_thread[bin][tid]+1.0;
	   //		printf("Intra beads are %d and %d and their distance is %lf\n",k,l,dist);
                        }
                    }
                }
                else {
                        if(dist<(rmax+0.0000001)&&(molid_bead[k]!=molid_bead[l])) {
                        bin=(int) (dist/dr);
                        gDeRinter_thread[bin][tid]=gDeRinter_thread[bin][tid]+1.0;
	//		printf("Inter beads are %d and %d and their distance is %lf\n",k,l,dist);
                        }
                }
                } //end of loop on pairs
            }  // end omp parallel
            
//	printf("Status 5  ok\n");
    /*Reduction*/
    for (idr=0;idr<nbin;idr++) {
        for (tid=0;tid<nthreads;tid++) {
            gDeRinter[idr] +=gDeRinter_thread[idr][tid];
            gDeRintra[idr] +=gDeRintra_thread[idr][tid];
       //     gDeRtot[idr] +=gDeRtot_thread[idr][tid];
                                        }
                                }
	
        iconf++;
    }

    printf("Analysis of %d configurations.\n",iconf);

    /*Normalisation of g(r) */
    // Normalisation of grbond
    norme=(float)(iconf*nmols*(nbeadspmol-1));
    for (idr=0;idr<nbin;idr++)
    gDeRbond[idr]=gDeRbond[idr]/norme;

    /*Normalisation of grbend */
    norme=(float)(iconf*nmols*(nbeadspmol-2));
    for (ii=0;ii<nbinctheta;ii++)
    gDeRbend[ii]=gDeRbend[ii]/norme;
    
    /*Normalisation of g(r) nonbonded */
//    for (idr=0;idr<nbin;idr++) {
//	}
    for (idr=0;idr<nbin;idr++) {
    vShell=4.0/3.0*acos(-1.)*(pow((idr+1)*dr,3)-pow(idr*dr,3));
    norme=(double)(iconf*npairs*(4.0/3.0)*acos(-1.)*(pow((idr+1)*dr,3)-pow(idr*dr,3))/(box[0]*box[1]*box[2]));
//	printf("Volume of shell is:%lf\n",vShell);
            gDeRinter[idr] = gDeRinter[idr]/norme;
            gDeRintra[idr] = gDeRintra[idr]/norme;
        //    gDeRtot[idr] /=norme;
    }  

 
	/* Impression */
	    sprintf(grbondName,"grbond_%d.out",z);
            grbondfile=fopen(grbondName,"w");
	if (grbondfile==NULL) 
		{
		printf("Can't open grbondfile \n");
			exit(1);
		} 
		sprintf(grbendName,"grbend_%d.out",z);
            grbendfile=fopen(grbendName,"w");
	if (grbendfile==NULL) 
		{
		printf("Can't open grbendfile \n");
			exit(1);
		} 
		sprintf(grintraName,"grintra_%d.out",z);
            grintrafile=fopen(grintraName,"w");
	if (grintrafile==NULL) 
		{
		printf("Can't open grintrafile \n");
			exit(1);
		} 
		sprintf(grinterName,"grinter_%d.out",z);
            grinterfile=fopen(grinterName,"w");
	if (grinterfile==NULL) 
		{
		printf("Can't open grinterfile \n");
			exit(1);
		} 
             printf("Writing  radial distribution functions\n");
//     for (idr=0;idr<(nbin-1);idr++) 
 //          fprintf(grtotfile,"%lf %lf\n",(idr-1)*dr+dr/2.0,gDeRtot[idr]);

            for (idr=0;idr<(nbin-1);idr++) 
            fprintf(grintrafile,"%lf %lf\n",(idr+0.5)*dr,gDeRintra[idr]);

            for (idr=0;idr<(nbin-1);idr++) 
            fprintf(grinterfile,"%d %lf %lf\n",idr,(idr+0.5)*dr,gDeRinter[idr]);
            
            for (idr=0;idr<(nbin-1);idr++) {
            ri=(idr+0.5)*dr;
            fprintf(grbondfile,"%lf %lf\n",ri,gDeRbond[idr]/(ri*ri));
            }

            for (ictheta=0;ictheta<(nbinctheta-1);ictheta++) 
            fprintf(grbendfile,"%lf %lf\n",(ictheta+0.5)*dctheta-1,gDeRbend[ictheta]);

   /* Free memory */
fclose(grinterfile);
fclose(grintrafile);
fclose(grbondfile);
fclose(grbendfile);
for (i=0;i<nbin;i++) {
//free(gDeRtot_thread[i]);
free(gDeRinter_thread[i]);
free(gDeRintra_thread[i]);
}
free(gDeRinter_thread);
free(gDeRintra_thread);
free(gDeRbond);
free(gDeRbend);
free(gDeRinter);
free(gDeRintra);
 return 0;
}
