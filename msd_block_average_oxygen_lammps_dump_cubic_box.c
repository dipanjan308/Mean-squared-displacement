#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define mol 2000 //total number of water molecules
#define N1 mol//number of atoms of type 1 (O)
#define N2 (2*mol)//number of atoms of type 2 (H)
#define datafile1 "pw_nvt1.dump" //path to the first datafile
#define datafile2 "pw_nvt2.dump" //path to the second datafile
#define max_time_f1 20000//maximum timestep for first data-file
#define max_time_f2 20000//maximum timestep for second data-file
#define max_time (max_time_f1+max_time_f2)//maximum timestep 
#define GAP 10//gap between two readings in data file
#define GAP_r 2//gap between two readings for calculation (actual gap is GAP*GAP_r)
#define max_array_size ((max_time/(GAP*GAP_r))*(N1+N2))//maximum array size required to store all frames
#define temperature 400//in K
#define density 1.3407//gm/cc
#define step_length 5//in fs

int timestep,atoms,atom_id,atom_type;
int index1[N1],index2[N2];
long mas=max_array_size;
double xlo,xhi,ylo,yhi,zlo,zhi,xt,yt,zt;

void get_particle_index_from_the_first_frame(const char *str)
{
	//char infile[100];
	char str1[100],str2[100],str3[100],str4[100],str5[100],str6[100],str7[100],str8[100],str9[100],str10[100];
	int i1,i2,i,ind,n;
	double x,y,z,vx,vy,vz;
	
	//sprintf(infile,"../../pw_nvt.dump");

	FILE *fp1;
	fp1=fopen(str,"r");
	
	i=0;
	{
		fscanf(fp1,"%s%s",str1,str2);
		fscanf(fp1,"%d",&timestep);
		fscanf(fp1,"%s%s%s%s",str1,str2,str3,str4);
		fscanf(fp1,"%d",&atoms);
		fscanf(fp1,"%s%s%s%s%s%s%s%s%s",str1,str2,str3,str4,str5,str6,str7,str8,str9);
		fscanf(fp1,"%lf%lf%lf",&xlo,&xhi,&xt);
		fscanf(fp1,"%lf%lf%lf",&ylo,&yhi,&yt);
		fscanf(fp1,"%lf%lf%lf",&zlo,&zhi,&zt);
		fscanf(fp1,"%s%s%s%s%s%s%s%s%s%s",str1,str2,str3,str4,str5,str6,str7,str8,str9,str10);
		i1=0;i2=0;
		ind=i/GAP;
		for(n=0;n<atoms;n++)
		{
			fscanf(fp1,"%d%d%lf%lf%lf%lf%lf%lf",&atom_id,&atom_type,&x,&y,&z,&vx,&vy,&vz);
			if(atom_type==1)
			{
				index1[i1]=atom_id-1;
				i1++;
			}
			else if(atom_type==2)
			{
				index2[i2]=atom_id-1;
				i2++;
			}
			else
			{
				printf("ERROR!\n");
				break;
			}
		}
   }
   fclose(fp1);
}

void import_position_file(const char *str, double *xcod, double *ycod, double *zcod, int start_frame, int maxt)
{
	//char infile[100];
	char str1[100],str2[100],str3[100],str4[100],str5[100],str6[100],str7[100],str8[100],str9[100],str10[100];
	int i,ind,ind2,n;
	double x,y,z,vx,vy,vz;
	//sprintf(infile,"../../pw_nvt.dump");
	FILE *fp1;
	fp1=fopen(str,"r");
	
	for(i=0;i<maxt;i=i+GAP)
	{
		fscanf(fp1,"%s%s",str1,str2);
		fscanf(fp1,"%d",&timestep);
		fscanf(fp1,"%s%s%s%s",str1,str2,str3,str4);
		fscanf(fp1,"%d",&atoms);
		fscanf(fp1,"%s%s%s%s%s%s%s%s%s",str1,str2,str3,str4,str5,str6,str7,str8,str9);
		fscanf(fp1,"%lf%lf%lf",&xlo,&xhi,&xt);
		fscanf(fp1,"%lf%lf%lf",&ylo,&yhi,&yt);
		fscanf(fp1,"%lf%lf%lf",&zlo,&zhi,&zt);
		fscanf(fp1,"%s%s%s%s%s%s%s%s%s%s",str1,str2,str3,str4,str5,str6,str7,str8,str9,str10);
		ind=i/(GAP*GAP_r);
		for(n=0;n<atoms;n++)
		{
			fscanf(fp1,"%d%d%lf%lf%lf%lf%lf%lf",&atom_id,&atom_type,&x,&y,&z,&vx,&vy,&vz);
			if(i%(GAP*GAP_r)==0)
			{
				ind2=atom_id-1+(start_frame+ind)*(N1+N2);
				xcod[ind2]=x;
				ycod[ind2]=y;
				zcod[ind2]=z;
			}
		}
   }
   fclose(fp1);
}
void unfold_trajectory(double *xcod, double *ycod, double *zcod)
{
        int i,j,ind,indp,t;
        t=floor(max_time/(GAP*GAP_r));
        //printf("%.8lf\t%.8lf\t%.8lf\n",box_x,box_y,box_z);
        for(i=1;i<t;i++)
        {
                for(j=0;j<N1+N2;j++)
                {
                        ind=j+i*(N1+N2);
                        indp=j+(i-1)*(N1+N2);
                        xcod[ind]=xcod[ind]-rint((xcod[ind]-xcod[indp])/box_x)*box_x;
                        ycod[ind]=ycod[ind]-rint((ycod[ind]-ycod[indp])/box_y)*box_y;
                        zcod[ind]=zcod[ind]-rint((zcod[ind]-zcod[indp])/box_z)*box_z;
                }
        }
}
double calculate_sliding_msd_each_particle(double *xcod, double *ycod, double *zcod, int block_size, int tau, int start_ind, int part_ind)
{
	int i,indi,indf;
	double delx,dely,delz,delsq=0.0;
	for(i=0;i<block_size/2;i++)
	{
		indi=start_ind+part_ind+(N1+N2)*i;
		indf=start_ind+part_ind+(N1+N2)*(i+tau);
		delx=xcod[indi]-xcod[indf];
		dely=ycod[indi]-ycod[indf];
		delz=zcod[indi]-zcod[indf];
		delsq+=(delx*delx+dely*dely+delz*delz);
	}
	return delsq;
}
void calculate_msd_sliding(double *xcod, double *ycod, double *zcod, int block_size, int tau, int start_ind, double delsq[N1])
{
	int i,ind;
	for(i=0;i<N1;i++)
	{
		ind=index1[i];
		delsq[i]=calculate_sliding_msd_each_particle(xcod,ycod,zcod,block_size,tau,start_ind,ind);
	}
}
double add_array(double arr[], int dim)
{
	int i;
	double sm=0.0;
	for(i=0;i<dim;i++)
	{
		sm+=arr[i];
	}
	return sm;
}
void calculate_msd_particle_average_all_tau(double *xcod, double *ycod, double *zcod, double msd_all_tau[], int block_size, int start_ind)
{
	int i;
	double delsq[N1],tmp;
	for(i=0;i<block_size/2;i++)
	{	
		calculate_msd_sliding(xcod,ycod,zcod,block_size,i,start_ind,delsq);
		tmp=add_array(delsq,N1);
		tmp=2.0*tmp/(block_size*N1);
		msd_all_tau[i]=tmp;
	}
}
void print_msd_data(double msd_data[], double msd_data2[], int block_size, int total_blocks)
{
	char outfile[100];
	int j;
	sprintf(outfile,"msd_slidingT%dden%3.4lfblock%dgap%d",temperature,density,block_size,GAP*GAP_r);
	FILE *fp;
	fp=fopen(outfile,"w");
	fprintf(fp,"#time(fs)\tMSD(A^2)\tError_MSD(A^2)\n");
	for(j=0;j<block_size/2;j++)
	{
		fprintf(fp,"%d\t%0.8e\t%0.8e\n",j*step_length*GAP_r,msd_data[j],sqrt(msd_data2[j]-msd_data[j]*msd_data[j])/sqrt(1.0*total_blocks));
	}
	fclose(fp);
}
void print_msd_data_each_block(double msd_data[], int block_size, int block_no)
{
	char outfile[100];
	int j;
	sprintf(outfile,"msd_slidingT%dden%3.4lfblock%dgap%d_%d",temperature,density,block_size,GAP*GAP_r,block_no);
	FILE *fp;
	fp=fopen(outfile,"w");
	fprintf(fp,"#time(fs)\tMSD(A^2)\n");
	for(j=0;j<block_size/2;j++)
	{
		fprintf(fp,"%d\t%0.8e\n",j*step_length*GAP_r,msd_data[j]);
	}
	fclose(fp);
}
void block_average(double *xcod, double *ycod, double *zcod, int block_size)
{
	int total_blocks,i,j,start_ind;
	double avg_msd[block_size/2],avg_msd2[block_size/2],msd_all_tau[block_size/2];
	total_blocks=(2*max_time)/(GAP*GAP_r*block_size);
	total_blocks--;
	printf("total_blocks=%d\n",total_blocks);
	for(i=0;i<block_size/2;i++)
	{
		avg_msd[i]=0.0;
		avg_msd2[i]=0.0;
	}
	for(i=0;i<total_blocks;i++)
	{
		start_ind=(N1+N2)*i*block_size/2;
		calculate_msd_particle_average_all_tau(xcod,ycod,zcod,msd_all_tau,block_size,start_ind);
		print_msd_data_each_block(msd_all_tau,block_size,i);
		for(j=0;j<block_size/2;j++)
		{
			avg_msd[j]+=msd_all_tau[j];
			avg_msd2[j]+=msd_all_tau[j]*msd_all_tau[j];
		}
	}
	for(j=0;j<block_size/2;j++)
	{
		avg_msd[j]=avg_msd[j]/(1.0*total_blocks);
		avg_msd2[j]=avg_msd2[j]/(1.0*total_blocks);
	}
	print_msd_data(avg_msd,avg_msd2,block_size,total_blocks);
}
void print_coordinates_vs_time(int part_ind, double *xcod, double *ycod, double *zcod)
{
	int i,ind;
	double x,y,z;
	for(i=0;i<2000;i++)
	{
		ind=index1[part_ind]+i*(N1+N2);
		x=xcod[ind];
		y=ycod[ind];
		z=zcod[ind];
		printf("%d\t%0.8e\t%0.8e\t%0.8e\n",i,x,y,z);
	}
}
int main()
{
	double *xcod,*ycod,*zcod;
	xcod = malloc(mas * sizeof(double));
	ycod = malloc(mas * sizeof(double));
	zcod = malloc(mas * sizeof(double));

	int ind_frame;
	get_particle_index_from_the_first_frame(datafile1);
	ind_frame=0;
	import_position_file(datafile1,xcod,ycod,zcod,ind_frame,max_time_f1);
	ind_frame=floor(max_time_f1/(GAP*GAP_r));
	import_position_file(datafile2,xcod,ycod,zcod,ind_frame,max_time_f2);
	//unfold_trajectory(xcod,ycod,zcod);
	int block_size;
	block_size=800;
	block_average(xcod,ycod,zcod,block_size);
	//print_coordinates_vs_time(1999,xcod,ycod,zcod);
	return (0);
}
