#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>


float starttime, endtime;
clock_t start, end;
double cpu_time_used;

void insertionSort(double arr[],int n,int p_8es[]);
void writeto(int npoints,int k_nn,int** knn,int size);
void readfrom(int npoints,int size,double** arxikos_pinakas,int diastaseis,int rank);
///piges=https://computing.llnl.gov/tutorials/mpi/

int main(int argc,char* argv[])
{
///gia na doulepsei sosta prepei npoints/size >= k_nn
///pisis na bgainei arkibis i dieresi npoints/size
    int rank,i,j,k,tag=0;
    int x,size;
    int file_free=0;
    signed int p;
    double key;
    int npoints=1000;
    int k_nn=10;
    int diastaseis=784;

    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    MPI_Request req1,req2;
    MPI_Request send_request1,send_request2;


///edw orizw enan 2-d pinaka poy exei tis 8eseis ka8e simiou dld 0-n/p , n/p-2n/p ktl
    int** pinakas_8esewn=malloc((npoints/size)*sizeof(int*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    int *b5=malloc((npoints/size)*(npoints/size)*sizeof(int));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints/size; i++)
    {
        pinakas_8esewn[i]=&b5[i*(npoints/size)];
    }

///arxikopoiisi pinaka_8esewn
    #pragma omp parallel for schedule(static)
    for(i=0; i<npoints/size; i++)
    {

        pinakas_8esewn[i][0]=(rank*(npoints/size));
        for(j=1; j<npoints/size; j++)
        {
            pinakas_8esewn[i][j]=pinakas_8esewn[i][j-1]+1;
        }
    }


///autos einai o pinakas poy exei tis 8eseis ka8e simiou alla 8a kanei rcv kiolas apo alles diergasies(idios me pinaka 8esewn arxika)
    int** pinakas_8esewn_sendrecv=malloc((npoints/size)*sizeof(int*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    int *b6=malloc((npoints/size)*(npoints/size)*sizeof(int));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints/size; i++)
    {
        pinakas_8esewn_sendrecv[i]=&b6[i*(npoints/size)];
    }

///autos einai o pinakas poy exei tis 8eseis ka8e simiou alla 8a kanei send kiolas se alles diergasies(idios me pinaka 8esewn arxika)
    int** pinakas_8esewn_send=malloc((npoints/size)*sizeof(int*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    int *b10=malloc((npoints/size)*(npoints/size)*sizeof(int));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints/size; i++)
    {
        pinakas_8esewn_send[i]=&b10[i*(npoints/size)];
    }

    ///arxikopiisi
    #pragma omp parallel for schedule(static)
    for(i=0; i<npoints/size; i++)
    {
        for(j=0; j<npoints/size; j++)
        {
            pinakas_8esewn_sendrecv[i][j]=pinakas_8esewn[i][j];
            pinakas_8esewn_send[i][j]=pinakas_8esewn[i][j];
        }
    }

///edw orizw enan 2-d pinaka gia ka8e diergasia poy 8a periexei ena kommati tou arxeiou
    double** arxikos_pinakas=malloc((npoints/size)*sizeof(double*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    double *b2=malloc((npoints/size)*diastaseis *sizeof(double));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints/size; i++)
    {
        arxikos_pinakas[i]=&b2[i*diastaseis];
    }

///edw orizw enan 2-d pinaka poy 8a exei ta dedomena poy ginontai rcv metaksi twn diergasiwn(arxika einai idios me ton arxiko_pinaka)
    double** send_recv=malloc((npoints/size)*sizeof(double*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    double *b4=malloc((npoints/size)*diastaseis *sizeof(double));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints/size; i++)
    {
        send_recv[i]=&b4[i*diastaseis];
    }

    ///edw orizw enan 2-d pinaka poy 8a exei ta dedomena poy stelnontai metaksi twn diergasiwn(arxika einai idios me ton arxiko_pinaka)
    double** send2=malloc((npoints/size)*sizeof(double*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    double *b9=malloc((npoints/size)*diastaseis *sizeof(double));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints/size; i++)
    {
        send2[i]=&b9[i*diastaseis];
    }
///kanw enan pinaka dist pou periexei tis apostaseis gia ka8e diergasia
    double** dist=malloc((npoints/size)*sizeof(double*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    double *b3=malloc((npoints/size)*(npoints/size)*sizeof(double));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints/size; i++)
    {
        dist[i]=&b3[i*(npoints/size)];
    }


///kanw enan pinaka dist2 pou periexei tis apostaseis gia ka8e diergasia me ta kainouria simia pou ir8an
    double** dist2=malloc((npoints/size)*sizeof(double*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    double *b8=malloc((npoints/size)*(npoints/size)*sizeof(double));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints/size; i++)
    {
        dist2[i]=&b8[i*(npoints/size)];
    }

    #pragma omp parallel
    {
///arxikopiisi tou dist2(epeidi einai sum=0)
        #pragma omp for schedule(static) nowait
        for(i=0; i<npoints/size; i++)
        {
            for(j=0; j<npoints/size; j++)
            {
                dist[i][j]=0;
            }
        }
///arxikopiisi tou dist2(epeidi einai sum=0)
        #pragma omp for schedule(static) nowait
        for(i=0; i<npoints/size; i++)
        {
            for(j=0; j<npoints/size; j++)
            {
                dist2[i][j]=0;
            }
        }

    }

///autos einai o pinakas me toys kontinoterous geitones
    int** knn=malloc((npoints/size)*sizeof(int*));;///xrisimopioume malloc(heap) gt perissoteros xwros apo to stack
    int *b7=malloc((npoints/size)*k_nn*sizeof(int));///kanoume malloc me ayuton ton tropo gia na exw ta dedomena se sinexomenes 8eseis stin mnimi
    for(i=0; i<npoints/size; i++)
    {
        knn[i]=&b7[i*(k_nn)];
    }
///telos "data segment"

///to diabasma tou arxeiou ginetai me blocking diadikasies
    if(size%2==0)
    {
        ///ka8e diergasia diabazei to antistoixo block kalontas tin readfrom me tin seira
        if(rank==0)
        {
            file_free=1;
        }
        else
        {
            MPI_Recv(&file_free,1,MPI_INT,rank-1,tag,MPI_COMM_WORLD,&status);
        }

        if (file_free==1) ///kalese tin sinartisi gia procces tou file
        {
            readfrom(npoints,size,arxikos_pinakas,diastaseis,rank);
        }
        ///dose adeia sto epomeno procces
        if(rank!=size-1)
        {
            MPI_Send(&file_free,1,MPI_INT,rank+1,tag,MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);


///oles oi diergasies: (analoga to rank tous exoyn kai sigekrimeno block(p.x rank=2 exei apo 2*n/p mexri 3*n/p)

        ///starttime = MPI_Wtick();///epeidi den bgazei tpt me auto
        start=clock();

///arxika o send_recv kai o send2 pinakas einai idios me ton arxiko
        #pragma omp parallel for schedule(static)
        for(i=0; i<(npoints/size); i++)
        {
            for(j=0; j<diastaseis; j++)
            {
                send_recv[i][j]=arxikos_pinakas[i][j];
                send2[i][j]=arxikos_pinakas[i][j];
            }
        }


///edw pernoyme tis apostaseis ka8e simiou apo to block simiwn poy exei ka8e diergasia

        for(k=0; k<npoints/size; k++)
        {
            for(i=0; i<(npoints/size); i++)
            {
                for(j=0; j<diastaseis; j++)
                {
                    if(k==i)
                    {
                        dist[k][i]=9999; ///gia na min bgoyn oloi oi geitones 0 kanoyme auto to trick
                        break;///kanoume break gia na ginei to i=1 kai na min baloume j fores tin idia timi sto dist[k][i]
                    }
                    else
                    {
                        dist[k][i]+=(arxikos_pinakas[k][j]-send_recv[i][j])*(arxikos_pinakas[k][j]-send_recv[i][j]);
                    }
                }
            }
        }


        for(i=0; i<(npoints/size); i++)
        {
            insertionSort(dist[i],(npoints/size),pinakas_8esewn[i]);
        }


        #pragma omp parallel for schedule(static)
        for(i=0; i<npoints/size; i++)
        {
            for(j=0; j<k_nn; j++)
            {
                knn[i][j]=pinakas_8esewn[i][j];
            }
        }

///telos protou iteration(se auto to simio exw tous 30 kontinoterous geitones twn prwtwn simiwn ka8e diergasias)

        int flag;

        for(x=0; x<size-1; x++)
        {
///topologia daktiliou

            if ( rank == 0 )
            {
                MPI_Irecv(*pinakas_8esewn_sendrecv,(npoints/size)*(npoints/size), MPI_INT,size - 1, tag+1, MPI_COMM_WORLD,&req1);
                MPI_Irecv(*send_recv,(npoints/size)*diastaseis, MPI_DOUBLE,size - 1, tag, MPI_COMM_WORLD,&req2);
                MPI_Isend(&send2[0][0],(npoints/size)*diastaseis, MPI_DOUBLE, (rank + 1) % size,tag,MPI_COMM_WORLD,&send_request1);
                MPI_Isend(&pinakas_8esewn_send[0][0],(npoints/size)*(npoints/size), MPI_INT, (rank + 1) % size,tag+1,MPI_COMM_WORLD,&send_request2);

                // MPI_Test(&send_request2, &flag,MPI_STATUS_IGNORE);///kanoume auta mesa sto while gt einai aneksartita apo tous bufer twn send kai rcv
                /*  while (!flag  && x==0)/// x=0 wste na ginei mia fora  gia ka8e diergasia kai oxi size-1
                  {

                      /* Do some work ...
                      MPI_Test(&send_request2,&flag,MPI_STATUS_IGNORE);
                  }*/
            }
            else if( rank != 0 )
            {
                MPI_Isend(&send2[0][0],(npoints/size)*diastaseis, MPI_DOUBLE, (rank + 1) % size,tag,MPI_COMM_WORLD,&send_request1);
                MPI_Irecv(*send_recv,(npoints/size)*diastaseis, MPI_DOUBLE,rank - 1, tag, MPI_COMM_WORLD,&req1);
                MPI_Irecv(*pinakas_8esewn_sendrecv,(npoints/size)*(npoints/size), MPI_INT,rank - 1, tag+1, MPI_COMM_WORLD,&req2);
                MPI_Isend(&pinakas_8esewn_send[0][0],(npoints/size)*(npoints/size), MPI_INT, (rank + 1) % size,tag+1,MPI_COMM_WORLD,&send_request2);

                /*   MPI_Test(&send_request2, &flag,MPI_STATUS_IGNORE);///kanoume auta mesa sto while gt einai aneksartita apo tous bufer twn send kai rcv
                   while (!flag && x==0)///x=0 wste na ginei mia fora kai oxi size-1
                   {

                       MPI_Test(&send_request2,&flag,MPI_STATUS_IGNORE);
                   }*/
            }


            MPI_Wait(&send_request1,MPI_STATUS_IGNORE);
            MPI_Wait(&send_request1,MPI_STATUS_IGNORE);
            MPI_Wait(&req1,MPI_STATUS_IGNORE);
            MPI_Wait(&req2,MPI_STATUS_IGNORE);


            ///efoson ginei  i antallagi oles oi diergasies antigrafoun ston send2 ton send_rcv wste na ginei sosta h epomenh antalagh
            #pragma omp parallel
            {


                #pragma omp for schedule(static)
                for(i=0; i<npoints/size; i++)
                {
                    for( j=0; j<diastaseis; j++)
                    {
                        send2[i][j]=send_recv[i][j];
                    }
                }

                ///omoia me akribws prin
                #pragma omp  for schedule(static)
                for(i=0; i<npoints/size; i++)
                {
                    for( j=0; j<npoints/size; j++)
                    {
                        pinakas_8esewn_send[i][j]=pinakas_8esewn_sendrecv[i][j];
                    }
                }

                #pragma omp for schedule(static)
                for(i=0; i<npoints/size; i++)
                {
                    for(j=0; j<npoints/size; j++)
                    {
                        dist2[i][j]=0;
                    }
                }
            }
            ///edw pernoyme tis apostaseis ka8e simiou apo to block simiwn poy exei ka8e diergasia efoson egine i antalagi

            for(k=0; k<npoints/size; k++)
            {
                for(int y=0; y<(npoints/size); y++)
                {
                    for(j=0; j<diastaseis; j++)///edw den xreazetai na poyme an k==y gt ston dist2 briskoume tin apostasi diaforetikwn simiwn akoma kai otan k==y
                    {
                        dist2[k][y]+=(arxikos_pinakas[k][j]-send_recv[y][j])*(arxikos_pinakas[k][j]-send_recv[y][j]);
                    }
                }
            }


            ///edw  elenxoume gia ka8e stixio tou dist2 an exei mikroteri apostasi apo ton dist poy exoume opote an isxiei auto, bazoume ton kainourio geitona, ston knn,
            ///stin 8esi tou paliou kai allazoume kai tin timi tou dist=dist2 etsi wste otan ginei o epomenos elenxos na ginei me tin updated timi tou dist kai kanoume kai sort

            for(k=0; k<npoints/size; k++)
            {
                for(int y=0; y<(npoints/size); y++)
                {
                    key=dist2[k][y];
                    for(p=k_nn-1; p>=0; p--) ///auto einai k_nn=30 giati o dist einai idi sorted ara elenxw mono tis protes 30 times tou dist
                    {
                        if(key<dist[k][p])
                        {
                            knn[k][p]=pinakas_8esewn_sendrecv[k][y];///edw kanoume update ton knn kai ton dist kai stin sinexeia kanoume sort ton knn me basi ton dist wste
                            dist[k][p]=key;                         /// na exoume tous k_nn mikroterous geitones se auksousa seira
                            insertionSort(dist[k],k_nn,knn[k]);
                            break;///break wste na mpainei mono mia fora o ka8e geitonas
                        }
                    }
                }
            }
        }
        ///se auto to simiw exw tous kontinoterous geitones olwn twn simwn stous knn twn diergasiewn(dld rank=0 knn gia 0 - n/p,rank=1 knn gia n/p - 2n/p ktl)

        end=clock();
        ///dinoume diadoxika adeia stis diergasies na grapsoun sto arxeio
        file_free=0;///to ksanakanoume 0 giati xrisimopii8ike prin gia to read
        if(rank==0)
        {
            file_free=1;
        }
        else
        {
            MPI_Recv(&file_free,1,MPI_INT,rank-1,tag,MPI_COMM_WORLD,&status);
        }

        if (file_free==1) ///kalese tin sinartisi gia procces tou file
        {
            writeto(npoints,k_nn,knn,size);
        }
        ///dose adeia sto epomeno procces
        if(rank!=size-1)
        {
            MPI_Send(&file_free,1,MPI_INT,rank+1,tag,MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        /// endtime = MPI_Wtick();


        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        ///double cpu_time_used2=(double)endtime-starttime;
        printf("rank %d took %lf seconds\n",rank,cpu_time_used);

        printf("\n");


    }
    else
    {
        printf("Oi diergasies prepei na einai artios ari8mos. Terminating.\n");
    }

    free(b2);
    free(b3);
    free(b4);
    free(b5);
    free(b6);
    free(b7);
    free(b8);
    free(b9);
    free(b10);
    free(arxikos_pinakas);
    free(send_recv);
    free(send2);
    free(pinakas_8esewn);
    free(pinakas_8esewn_sendrecv);
    free(pinakas_8esewn_send);
    free(dist2);
    free(dist);
    free(knn);

    MPI_Finalize();


    return 0;
}

///se autin tin sinartisi kanoume sort ton array(dist)(me auksousa diataksi) kai ton pinaka_8esewn me basi ton array
///kai stin sinexeia xrisimopieite gia na kanei sort ton knn me basi ta prota k_nn stixia tou array(dist)
void insertionSort(double arr[],int n,int p_8es[])
{
    int i,j,key2;
    double key;

    for (i = 1; i < n; i++)
    {
        key = arr[i];
        key2 = p_8es[i];
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

///grafoume se arxeio ta apotelesmata
void writeto(int npoints,int k_nn,int **knn,int size)
{


    FILE* myfile;
    int i,j;

    myfile=fopen("data2.txt","a");
    if (!myfile)
    {
        printf("Unable to open file for write!");
        exit(1);
    }
    //char buffer[npoints/size][k_nn];

    for(i=0; i<npoints/size; i++)
    {
        for(j=0; j<k_nn; j++)
        {
            //buffer[i][j]=knn[i][j] + '0';
            //fwrite(*buffer,sizeof(char),(npoints/size)*k_nn,myfile);///grafw ton knn sto arxeio
            //fprintf(myfile,"%c ",buffer[i][j]);
            fprintf(myfile,"%d ",knn[i][j]);
        }
        fprintf(myfile,"\n");
    }

    fclose(myfile);

}

///diabazoume apo to arxeio gia na kataneimoume se ka8e diergasia ena kommati
void readfrom(int npoints,int size,double** arxikos_pinakas,int diastaseis,int rank)
{

    FILE* myfile;
    int i,j;

    myfile=fopen("data3.bin","rb");
    if (!myfile)
    {
        printf("Unable to open file sto read!");
        exit(1);
    }

    fseek(myfile,rank*(npoints/size)*diastaseis*sizeof(double),SEEK_SET);
    fread(*arxikos_pinakas,sizeof(double),(npoints/size)*diastaseis,myfile);///diabazw olo ton pinaka

    fclose(myfile);

}


