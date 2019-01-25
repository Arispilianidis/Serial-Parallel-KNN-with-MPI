//#include "mpi.h" siriakooooooo
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

struct timeval startwtime, endwtime;
double seq_time;

void insertionSort(double arr[],int n,int p_8es[]);

///piges=https://computing.llnl.gov/tutorials/mpi/…

int main(int argc,char* argv[])
{
    gettimeofday( &startwtime, NULL );
///gia na doulepsei sosta prepei npoints/size >=k_nn SOSOSSOSSOSOSOSOSOOSOSOSOSOSOSSOOSSOSOS
    int rank,i,j,k,k_nn=10,tag=0;
    int size;
    double key;
    int npoints=1000;
    int diastaseis=20;
///omp_set_num_threads(omp_get_num_procs());///4 threads osoi oi dia8esimoi epeksergastes

    /**printf("posa simia 8es?\n");
    scanf(" %d",&npoints);
    printf("poses diastaseis exei to ka8e simio?\n");
    scanf(" %d",&diastaseis);
    printf("posous kontinoterous geitones psaxneis?\n");
    scanf(" %d",&k_nn);*/

///edw orizw enan 2-d pinaka poy exei tis 8eseis ka8e simiou dld 0-n/p , n/p-2n/p ktl
    int** pinakas_8esewn=malloc(npoints*sizeof(int*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    int *b5=malloc(npoints*npoints*sizeof(int));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints; i++)
    {
        pinakas_8esewn[i]=&b5[i*(npoints)];
    }

///arxikopoiisi pinaka_8esewn
    for(i=0; i<npoints; i++)
    {
        pinakas_8esewn[i][0]=0;

///#pragma omp parallel for
        for(j=1; j<npoints; j++)
        {
            pinakas_8esewn[i][j]=pinakas_8esewn[i][j-1]+1;
        }
    }

///kanw enan pinaka dist pou periexei tis apostaseis gia ka8e diergasia
    double** dist=malloc(npoints*sizeof(double*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    double *b3=malloc(npoints*npoints*sizeof(double));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints; i++)
    {
        dist[i]=&b3[i*(npoints)];
    }

///arxikopiisi tou distance(epeidi einai sum=0)
///#pragma omp parallel for
    for(i=0; i<npoints; i++)
    {
        for(j=0; j<npoints; j++)
        {
            dist[i][j]=0;
        }
    }

///autos einai o pinakas me toys kontinoterous geitones(k_nn gt ena simio einai gia auto pou psaxnoyme tous geitones)
    int** knn=malloc(npoints*sizeof(int*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    int *b7=malloc(npoints*k_nn*sizeof(int));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints; i++)
    {
        knn[i]=&b7[i*(k_nn)];
    }

///sto rank=0 diabazw olo to arxeio kai diamirazei ta kommatia sta ypolipa gt den katafera na to kanw me fseek() prin na kalesw tin fread.

    FILE* myfile;

    double** pinakas=malloc((npoints)*sizeof(double*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    double *b=malloc((npoints)*diastaseis *sizeof(double));///kanoume malloc me ayuton ton tropo gia na exw sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints; i++)
    {
        pinakas[i]=&b[i*diastaseis];
    }

    myfile=fopen("data.bin","rb");
    if (!myfile)
    {
        printf("Unable to open file!");
        return 1;
    }

    fread(*pinakas,sizeof(double),npoints*diastaseis,myfile);///diabazw olo ton pinaka
    fclose(myfile);

////////////////////////////////////////////////////starttime = MPI_Wtick();//auta na ginooun opws sto openmp gt den exoume mpi edw

    for(k=0; k<npoints; k++)
    {
        for(i=0; i<npoints; i++)
        {
            for(j=0; j<diastaseis; j++)
            {
///gia na min bgoyn oloi oi geitones 0 kanoyme auto to trick
                if(k==i)
                {
                    dist[k][i]=9999;
                    break;///kanoume break gia na ginei to i=1 kai na min baloume j fores tin idia timi sto dist[k][i]
                }
                else
                {
                  //  if (dist[k][i]==0) den einai tetragonikos-simmetrikos gia na to kanw auto
                    //{
                        dist[k][i]+=(pinakas[k][j]-pinakas[i][j])*(pinakas[k][j]-pinakas[i][j]);

                     /*   if (dist[i][k]==0)
                        {
                            dist[i][k]=dist[k][i];
                        }
                       }
                        else{

                        break;
                    }*/
                }
            }
        }
    }

///gt den ekana  qsort? -> gt i8ela na allazw kai to pinakas 8esewn analoga to dist apo oti 8imamai
for(i=0; i<npoints; i++)
{
    insertionSort(dist[i],npoints,pinakas_8esewn[i]);
}

for(i=0; i<npoints; i++)
{
    for(j=0; j<k_nn; j++)
    {
        knn[i][j]=pinakas_8esewn[i][j];
    }
}

free(b);
free(pinakas);

gettimeofday( &endwtime, NULL );
seq_time = (double)( ( endwtime.tv_usec - startwtime.tv_usec ) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec );

printf("THE SERIAL took %.2f seconds\n",seq_time);
printf("\n");


for(i=0; i<10; i++)
{
    printf("knn for %d:",i);
    for(j=0; j<k_nn; j++)
    {
        printf(" %d ",knn[i][j]);
    }
    printf("\n");
}

free(b3);
free(b5);
free(b7);

free(pinakas_8esewn);
free(dist);
free(knn);



return 0;
}

//////////////////////////////////////////////////////////////////////////////endtime = MPI_Wtick();

void insertionSort(double arr[],int n,int p_8es[])
{
    int i,j;
    double key,key2;

    for (i = 1; i < n; i++)
    {
        key = arr[i];
        key2=p_8es[i];
        j = i-1;

        /* Move elements of arr[0..i-1], that are
        greater than key, to one position ahead
        of their current position*/
        while (j >= 0 && arr[j] > key)
        {
            arr[j+1] = arr[j];
            p_8es[j+1]=p_8es[j];
            j = j-1;
        }
        arr[j+1] = key;
        p_8es[j+1]=key2;
    }
}
