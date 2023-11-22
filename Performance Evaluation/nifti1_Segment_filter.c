#define _GNU_SOURCE
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <assert.h>
#include <dirent.h>
#include "nifti1_read_write.h"


/*
Floating-point Max-tree algorithm

Michael Wilkinson

*/
#define BOTTOM (-1)

#define PI 3.14159265358979323846
#define false 0
#define true  1
#define MAXTHREADS 128
#define MAXPATHLEN 500
int MULFACTOR = 1;     

#define CONNECTIVITY  6
#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))
#define LWB(self) (size2D*(((self)*depth)/nthreads))
#define UPB(self) (size2D*(((self+1)*depth)/nthreads))
typedef short bool;
typedef unsigned char ubyte;

int ParCount=0;
int nthreads;
pthread_t threadID[MAXTHREADS];

int width, height, depth, size;  /* precondition: width <= size/nthreads */
int size2D;
double lambda;
clock_t start,built, merged,finish;
struct tms tstruct;

bool weightedfilter = true;

typedef MY_DATATYPE Pixel;

Pixel *gval=NULL, *out=NULL;

typedef struct MaxNode  
{ 
  int parent;
  int TumourLabel;
  bool filtered; /* indicates whether or not the filtered value is OK */
  Pixel outval;
} MaxNode;

#define bottom (-1)


MaxNode *node; 

/************************** safe malloc and calloc ************************************/

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void *SafeMalloc(int n) {
  void *ptr;
  pthread_mutex_lock(&mutex);
  ptr = malloc(n);
  if (ptr==NULL) {
    fprintf (stderr, "Error: out of memory. "
	     "Could not allocate %d bytes\n", n);
  }
  pthread_mutex_unlock(&mutex);
  return ptr;
}

void *SafeCalloc(int nmemb, int size) {
  void *ptr;
  pthread_mutex_lock(&mutex);
  ptr = calloc(nmemb, size);
  if (ptr==NULL) {
    fprintf (stderr, "Error: out of memory. Could not "
	     "allocate %d bytes\n", nmemb*size);
  }
  pthread_mutex_unlock(&mutex);
  return ptr;
}

/************************** semaphore ************************************/

pthread_mutex_t samut[MAXTHREADS];
pthread_cond_t  sacv[MAXTHREADS];
int             saval[MAXTHREADS];

void Psa(int p) {
  pthread_mutex_lock(&samut[p]);
  while (saval[p] <= 0)
    pthread_cond_wait(&sacv[p], &samut[p]);
  saval[p]--;
  pthread_mutex_unlock(&samut[p]);
}

void Vsa(int p) {
  pthread_mutex_lock(&samut[p]);
  saval[p]++;
  pthread_mutex_unlock(&samut[p]);
  pthread_cond_broadcast(&sacv[p]);
}

/************************** barrier ************************************/

pthread_mutex_t barriermutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t barriercv = PTHREAD_COND_INITIALIZER;
int barcnt = 0;

void Barrier(int self) {
  pthread_mutex_lock(&barriermutex);
  barcnt++;
  if (barcnt == nthreads) {
    barcnt = 0;  /* for reuse of routine */
    pthread_cond_broadcast (&barriercv);  
  } else {
    pthread_cond_wait (&barriercv, &barriermutex);
  }
  pthread_mutex_unlock(&barriermutex);
}

/************************** level root ************************************/

ulong levroot(ulong x) {
  int r=x, y ,gv=gval[x];
  if (r==bottom) return bottom;
  while ((node[r].parent!=bottom) && (gv==gval[node[r].parent]))
      r = node[r].parent;

  while (x!=r){
    y=node[x].parent;
    node[x].parent=r;
    x=y;
  }
  return r;
}

ulong Par(ulong x) {
  ParCount++;
  return levroot(node[x].parent);
}

void levrootfix(ulong lwb, ulong upb) {
  ulong z, u, x;

  for (x=lwb; x<upb;x++) {
    u = levroot(x);
    if (x!=u) node[x].parent=u;
    else {
      z = Par(x); 
      node[x].parent=z;
    }
  }
}

typedef struct{
  int  size, maxsize;
  int *stack;  
} pStack;


pStack *pStackCreate(long maxsize){
  pStack *newStack = (pStack *) malloc(sizeof(pStack));
  newStack->size = 0;
  newStack->stack = (int *)malloc((maxsize)*sizeof(int));
  newStack->maxsize=maxsize;
  return newStack;
}

void pStackDelete(pStack *oldstack){
  free(oldstack->stack);
  free(oldstack);
}

#define IsEmpty(stack)        ((stack->size)==0)
#define IsFull(stack)         (stack->size == stack->maxsize)
#define pStackTop(stack)       (stack->stack[stack->size-1])
#define pStackPop(stack)       (stack->stack[--stack->size])
#define pStackPush(stack,elem)      (stack->stack[stack->size++]=elem)

typedef struct{
  int size,maxsize;
  int *queue;
} pQueue;

pQueue *pQueueCreate(long maxsize){
  pQueue *newQueue = (pQueue *) malloc(sizeof(pQueue));
  newQueue->size = 0;
  newQueue->queue = (int *)malloc((maxsize+1)*sizeof(int));
  newQueue->maxsize=maxsize;
  return newQueue;
}


#define pQueueFront(queue)       (queue->queue[1])

void pQueueDelete(pQueue *oldqueue){
  free(oldqueue->queue);
  free(oldqueue);

}

int pQueuePop(pQueue *queue, Pixel *priority){
  int outval = queue->queue[1];
  int current = 1,
    moved;
  Pixel curval;

  moved = queue->queue[queue->size];
  queue->size--;
  curval = priority[moved];

  while ( ((current*2<=queue->size) &&
	   (curval< priority[queue->queue[current*2]]))
	  ||
	  ((current*2+1<=queue->size) &&
	   (curval< priority[queue->queue[current*2+1]]))
	  ){
    if ((current*2+1<=queue->size) && 
	(priority[queue->queue[current*2]]< 
	 priority[queue->queue[current*2+1]])){
      queue->queue[current]= queue->queue[current*2+1];
      current+=current+1;
    } else {
      queue->queue[current]= queue->queue[current*2];
      current+=current;
    }
  }
  queue->queue[current]=moved;
  return outval;
}

void pQueuePush(pQueue *queue, Pixel *priority, int pixpos){
  long current;
  Pixel curval = priority[pixpos];
  queue->size++;
  current=queue->size;
  
  while ((current/2 !=0) && (priority[queue->queue[current/2]]<curval)){
    queue->queue[current]= queue->queue[current/2];
    current=current/2;
  }

  queue->queue[current]=pixpos;
}

int GetNeighbors(int p, int x, int y, int z, 
		 int *neighbors, int width, int height, int depth,
		 int size2D, int lwb, int upb)
{

   int n=0;


   if (x<(width-1))       neighbors[n++] = p+1;
   if (y>0)               neighbors[n++] = p-width;
   if (x>0)               neighbors[n++] = p-1;
   if (y<(height-1))      neighbors[n++] = p+width;
   if (depth>1) {
      if (z>lwb)          neighbors[n++] = p-size2D;
      if (z<(upb-1))           neighbors[n++] = p+size2D;
   }
   return(n);
} /* GetNeighbors */


void Flood(int self, pQueue *queue, pStack *stack, 
	   Pixel *gval, int width, int height, int depth,
	   MaxNode *node){
  int neighbors[CONNECTIVITY];
  int size = width*height*depth;
  int size2D = width*height;
  int lwb, upb;

  int xm,i, x,y,z, nextpix, p, q, numneighbors, oldtop;

  lwb = LWB(self); upb = UPB(self);
  xm = lwb;
  for (x=xm; x<upb; x++) { 
    node[x].parent = bottom;
    node[x].filtered = false;
    if (gval[xm]>gval[x]) xm = x;
  }

  pQueuePush(queue,gval,xm);
  nextpix = xm;

  node[nextpix].parent = nextpix;
  pStackPush(stack,nextpix);
  node[nextpix].TumourLabel = 0;
  

  do{ 
    p=nextpix;
      
    x = p % width; 
    y = (p % size2D) / width;
    z = p / size2D;

    numneighbors = GetNeighbors(p, x, y, z, neighbors, width, height,
				depth, size2D, lwb/size2D, upb/size2D);    
    for (i=0; i<numneighbors; i++) {
      q = neighbors[i];
      if (node[q].parent == BOTTOM){
	node[q].parent=q;
	node[q].TumourLabel=0;
	pQueuePush(queue,gval,q);
        if (gval[q]>gval[p]){
	  break;
	}
      }
    }
    nextpix = pQueueFront(queue);

    if (gval[nextpix]>gval[p]){
      /* start processing at a higher level */      
      pStackPush(stack,nextpix);
    } else { /* nextpix == p */
      p=pQueuePop(queue,gval); /* does not change p but removes from queue */
      if (p!=pStackTop(stack)){
	node[p].parent = pStackTop(stack);
	/*node[pStackTop(stack)].TumourLabel++;*/
      }

      nextpix = pQueueFront(queue); /* new front of queue */
      /* if queue is empty, nextpix=p */

      if (gval[nextpix]<gval[p]){
	/* descend to lower level */	
	do { 
	  oldtop=pStackPop(stack);
	  
	  if (gval[pStackTop(stack)]<=gval[nextpix]) 
	    break;
	  node[oldtop].parent = pStackTop(stack);
	  node[pStackTop(stack)].TumourLabel += node[oldtop].TumourLabel;
	} while (1);
	
	
	if (gval[pStackTop(stack)]< gval[nextpix]){
	  pStackPush(stack,nextpix);
	}
 
	node[oldtop].parent = pStackTop(stack);
	node[pStackTop(stack)].TumourLabel += node[oldtop].TumourLabel;
	/* removed everything from stack down to gval[nextpix]*/
      }
    }
  } while (!IsEmpty(queue));

  /* fix remainder of stack */

  oldtop=pStackPop(stack);
  while (!IsEmpty(stack)) { 
    node[oldtop].parent = pStackTop(stack);
    node[pStackTop(stack)].TumourLabel += node[oldtop].TumourLabel;
    oldtop=pStackPop(stack);
  }
    
  node[xm].parent=bottom;
}

void SetTumourLabels(MaxNode *node, char *CSVfile)
{
  int i,total=0;
  FILE *fp;
  size_t linelength=0; 
  char *linebuffer=  malloc(linelength*sizeof(char));
  int node_id,parent_id,label;
  ssize_t charsfound;  
  
  fp = fopen(CSVfile,"r");
  if (fp == NULL) {
        fprintf(stderr, "\nError opening header file %s\n",CSVfile);
        exit(1);	
  }
    
  charsfound = getline(&linebuffer, &linelength, fp);
  while ((charsfound = getline(&linebuffer, &linelength, fp)) != -1){
    sscanf(linebuffer,"%d,%d,%d",&node_id,&parent_id,&label);

    node[node_id].TumourLabel=label;
    total+=label;
  }
  free(linebuffer);
  printf("Total labels = %d\n",total);
}

void MaxTreeTumourLabelFilter(int self, double lambda)
{
   ulong v, u, w, parent, lwb, upb;
   Pixel val;
 
   lwb = LWB(self); upb = UPB(self);
   for (v=lwb; v<upb; v++) {
      if (!node[v].filtered) {
         w = v;
         parent = node[w].parent;
         while ((parent != bottom) && (!node[w].filtered) && 
            ((gval[w] == gval[parent]) || (node[w].TumourLabel == 0))) {
                  w = parent;
                  parent = node[w].parent;
         }
         if (node[w].filtered) val = node[w].outval;
         else if (node[w].TumourLabel) val = 1; 
         else val = 0; /* criterion cannot be satisfied */
         u = v;
         while (u!=w) {
           if ((lwb<=u) && (u<upb)){ 
	     node[u].outval = val;
	     node[u].filtered=true;
	   }

           u = node[u].parent;
           }
         if ((lwb<=w) && (w<upb)){ 
	   node[w].outval = val;
	   node[w].filtered=true;
	 }
       }
      out[v]=node[v].outval;
   }
} /* MaxTreeTumourLabelFilter */

void MaxTreeTumourLabelFilterWeighted(int self, double lambda, pStack *stack)
{
   ulong v, u, w, parent, lwb, upb;
   Pixel val;
   
   stack->size=0;
   
   lwb = LWB(self); upb = UPB(self);
   for (v=lwb; v<upb; v++) {
      if (!node[v].filtered) {
         w = v;
         parent = node[w].parent;
         while ((parent != bottom) && (!node[w].filtered)) {
	   pStackPush(stack,w);
	   w = parent;
	   parent = node[w].parent;
         }
         if (node[w].filtered)
	   val = node[w].outval;
         else {
	   node[w].outval=0; /* w is root*/
	   val = 0;
	   node[w].filtered = true;
	 }
	 while (!IsEmpty(stack)){
	   w=pStackPop(stack);	   
	   val += node[w].TumourLabel;
	   if ((lwb<=w) && (w<upb)){ 
	     node[w].outval = val;
	     node[w].filtered=true;
	   }
	 }
       }
      out[v]=node[v].outval;
   }
} /* MaxTreeTumourLabelFilter */





void Connect(ulong x, ulong y) {

  ulong area = 0, area1 = 0;
  ulong h, z;

  x = levroot(x);
  y = levroot(y);
  if (gval[x] < gval[y]) {
    h=x; x=y; y=h;
  }
  while ((x!=y) && (y!=bottom)) {
    z = Par(x);
    if ((z!=bottom) && (gval[z]>=gval[y])) {
      node[x].TumourLabel += area;
      x = z;
    } else {
      area1 = node[x].TumourLabel + area;
      area = node[x].TumourLabel ;
      node[x].TumourLabel =area1;
      node[x].parent = y;
      x = y;
      y = z;
    }
  }
  if (y==bottom) {
    while(x!=bottom) {
      node[x].TumourLabel += area;
      x = Par(x);
    }
  }
}

void Fuse2(int self, int i) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
{
  ulong lwb, upb;
  ulong p, xm, x, y;
   
  lwb = LWB(self);
  upb = UPB(self+2*i-1); 
  if (upb>size) upb=size;

  /* get the horizontal boundary */
  xm = LWB(self+i); 

  x = xm;
  p = x % width;
  if ((p>0) && (x-1>=lwb)) {
     y = x-1;
     Connect(x, y);
  }  

  for (p=0; p<width && x<upb; p++){
      if (x>=lwb+width) {
         y= x-width;
         Connect(x, y);
      }
      x++;
  }

  if (depth>1) {
    x = xm;
    for (p=0; p<size2D && x<upb; p++){
      if (x>=lwb+size2D) {
	y= x-size2D;
	Connect(x, y);
      }
      x++;
    }
  }

  /* levrootfix(lwb, upb); */
}


void Fuse(int self, int i) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
{
  ulong p, q, x, y, z;

  Pixel prevmin, curmin, nextmin;

  /* get the horizontal boundary */ 

  p  = LWB(self+i);
  q = p - size2D;

  x = p % width;
  y = (p /width) % height;
  z = p /size2D;
  
  /*  printf("Region %d merger with %d: (%d,%d,%d)\n",self, self+i, x,y,z);*/

  for ( y = 0 ; y < height ; y++ ){
    Connect(p, q);
    prevmin = MIN(gval[p],gval[q]);
    p++;
    q++;
    curmin = MIN(gval[p],gval[q]);
    for ( x = 1 ; x < width-1 ; x++, p++, q++ ){
      nextmin = MIN(gval[p+1],gval[q+1]);
      if ((curmin > prevmin) && (curmin>=nextmin)){
	Connect(p, q);
      }
      prevmin = curmin;    
      curmin = nextmin;
    }  
    p++;
    q++;
    Connect(p, q);

  }

  /* levrootfix(lwb, upb); */
}


/**************** Concurrent construction and filter of Maxtree  ***********/

typedef struct { 
  int self;
  pQueue *thisqueue;
  pStack *thisstack;
  char *CSVfname;
} ThreadData;

ThreadData *MakeThreadData(int numthreads, char *CSVfname){
  ThreadData *data = malloc(numthreads *sizeof(ThreadData));
  int i;

  for (i=0; i<numthreads; i++){
    data[i].self=i;
    data[i].thisstack= pStackCreate( UPB(i)-LWB(i));
    data[i].thisqueue= pQueueCreate( UPB(i)-LWB(i));
    data[i].CSVfname= CSVfname;
    printf("Thread %d, LWB = %d, UPB = %d\n",i,LWB(i),UPB(i));
  }    
  return(data);
}

void FreeThreadData(ThreadData *data, int numthreads)
{
  int i;
  for (i=0;i<numthreads; i++){
    pQueueDelete(data[i].thisqueue);
    pStackDelete(data[i].thisstack);
  }
  free(data);
}


void *ccaf(void *arg) {
  ThreadData *thdata = (ThreadData *) arg;
  int self = thdata->self, q, i;
  ulong x, area=0;


  Flood( self, thdata->thisqueue, thdata->thisstack, 
	 gval, width, height, depth, node);


  i = 1;
  q = self;
  if (self==0)
    built = times(&tstruct);

  while ((self+i<nthreads) && (q%2 == 0)) {
    Psa(self+i);  /* wait to glue with righthand neighbor */
    Fuse(self, i);
    i = 2*i;
    q = q/2;
  }
  if (self != 0) {
    Vsa(self);  /* signal lefthand neighbor */
  }
  Barrier(self);
  if (self==0){
    merged = times(&tstruct);
    SetTumourLabels(node, thdata->CSVfname);
  }
  if (weightedfilter)
    MaxTreeTumourLabelFilterWeighted(self, lambda, thdata->thisstack );
  else
    MaxTreeTumourLabelFilter(self, lambda);
    
  return NULL;
}


void BuildTreeAndFilter(ThreadData *thdata, int nthreads) {
  int thread;
  for (thread=0; thread<nthreads; ++thread) {
    pthread_create(threadID+thread, NULL, ccaf, (void *) (thdata + thread));
  }
  for (thread=0; thread<nthreads; ++thread) {
    pthread_join(threadID[thread], NULL);
  }
}

int main (int argc, char *argv[]) {

  char *in_dir_path, *out_dir_path, *imgfname, *CSVfname, *outfname = "out.nii";
  ThreadData *thdata;
  nifti_1_header header;
  int r;
  ulong i;
  long tickspersec = sysconf(_SC_CLK_TCK);  
  float musec;

  if (argc<4)
    {
      printf("Usage: %s <input directory> <output directory> <nthreads>\n", argv[0]);
      exit(0);
   }

   in_dir_path = argv[1];
   out_dir_path = argv[2];
   nthreads = MIN(atoi(argv[3]), MAXTHREADS);
   
   DIR *in_dir;
   struct dirent *entry;
   in_dir = opendir(in_dir_path);
   
   if (in_dir == NULL) {
      perror("opendir");
      exit(EXIT_FAILURE);
   }
   
   while ((entry = readdir(in_dir)) != NULL) {
      if (strstr(entry->d_name, ".nii") == NULL) {
         continue;
      }

      char nii_path[MAXPATHLEN];
      sprintf(nii_path, "%s/%s", in_dir_path, entry->d_name);
      
      char base_name[MAXPATHLEN];
      strcpy(base_name, entry->d_name);
      base_name[strlen(base_name)-4] = '\0';
      
      char csv_path[MAXPATHLEN];
      char out_path[MAXPATHLEN];
      sprintf(csv_path, "%s/%s.csv", in_dir_path, base_name);
      sprintf(out_path, "%s/%s.nii", out_dir_path, base_name);
      if (read_nifti_file(nii_path, nii_path, &header, &gval)) {
         fprintf(stderr, "Error reading input file: %s\n", nii_path);
         continue;
      }

      width=header.dim[1];
      height=header.dim[2];
      depth=header.dim[3];
      size = width*height*depth;
      size2D = width*height;

      printf("Filtering image '%s' using attribute area with lambda=%f\n", nii_path, lambda);
      printf("Image: Width=%d Height=%d Depth=%d\n", width, height, depth);
      printf ("nthreads: %d\n", nthreads);

      node = calloc((size_t)size, sizeof(MaxNode));
      if (node==NULL) {
         fprintf(stderr, "out of memory! \n");
         free(gval);
         continue;
      }

      out =  malloc(size*sizeof(Pixel));
      if (out==NULL) {
        fprintf(stderr, "Can't create output image! \n");
        free(node);
        free(gval);
        continue;
      }

      thdata = MakeThreadData(nthreads, csv_path); 
      printf("Data read, start filtering.\n");
      start = times(&tstruct);

      BuildTreeAndFilter(thdata,nthreads);

      finish = times(&tstruct);
      musec = (float)(built - start)/((float)tickspersec);

      printf("build time: %f s\n",musec);
      musec = (float)(merged - built)/((float)tickspersec);

      printf("merge time: %f s\n",musec);
      musec = (float)(finish - merged)/((float)tickspersec);

      printf("filtering time: %f s\n",musec);
      musec = (float)(finish - start)/((float)tickspersec);

      printf("wall-clock time: %f s\n",musec);
   

      printf("Parcount = %d\n",ParCount); 

     if (write_nifti_file(out_path, out_path,&header, out)) {
      fprintf(stderr, "Error writing output file: %s\n", out_path);
      continue;}
      free(out);
      if (r)  printf("Filtered image written to '%s'\n", outfname);
   
      FreeThreadData(thdata,nthreads);
      free(node);
      free(gval);
     
   
}
closedir(in_dir);
return 0;
} /* main */