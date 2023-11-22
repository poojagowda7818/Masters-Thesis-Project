/*********************************************************************
 *
 * Very simple code snippets to read/write nifti1 files
 * This code is placed in the public domain.
 *
 * If you are the type who doesn't want to use a file format unless
 * you can write your own i/o code in less than 30minutes, this
 * example is for you.
 *
 * This code does not deal with wrong-endian data, compressed data,
 * the new qform/sform orientation codes, parsing filenames, volume-
 * wise or timecourse-wise data access or any of a million other very useful
 * things that are in the niftilib i/o reference libraries.
 * We encourage people to use the niftilib reference library and send
 * feedback/suggestions, see http://niftilib.sourceforge.net/
 * But, if that is too much to tackle and you just want to jump in, this
 * code is a starting point.
 * This code was written for maximum readability, not for the greatest
 * coding style.
 *
 *
 * If you are already a little familiar with reading/writing Analyze
 * files of some flavor, and maybe even have some of your own code, here
 * are the most important things to be aware of in transitioning to nifti1:
 *
 * 1. nii vs .hdr/.img
 *      nifti1 datasets can be stored either in .hdr/.img pairs of files
 *      or in 1 .nii file.  In a .nii file the data will start at the byte
 *      specified by the vox_offset field, which will be 352 if no extensions
 *      have been added.  And, nifti1 really does like that magic field set
 *      to "n+1" for .nii and "ni1" for .img/.hdr
 *
 * 2. scaling
 *      nifti1 datasets can contain a scaling factor.  You need to check the
 *      scl_slope field and if that isn't 0, scale your data by
 *      Y * scl_slope  + scl_inter
 *
 * 3. extensions
 *      nifti1 datasets can have some "extension data" stuffed after the
 *      regular header.  You can just ignore it, but, be aware that a
 *      .hdr file may be longer than 348 bytes, and, in a .nii file
 *      you can't just jump to byte 352, you need to use the vox_offset
 *      field to get the start of the image data.
 *
 * 4. new datatypes
 *      nifti1 added a few new datatypes that were not in the Analyze 7.5
 *      format from which nifti1 is derived.  If you're just working with
 *      your own data this is not an issue but if you get a foreign nifti1
 *      file, be aware of exotic datatypes like DT_COMPLEX256 and mundane
 *      things like DT_UINT16.
 *
 * 5. other stuff
 *     nifti1 really does like the dim[0] field set to the number of
 *     dimensions of the dataset.  Other Analyze flavors might not
 *     have been so scrupulous about that.
 *     nifti1 has a bunch of other new fields such as intent codes,
 *     qform/sform, etc. but, if you just want to get your hands on
 *     the data blob you can ignore these.  Example use of these fields
 *     is in the niftilib reference libraries.
 *
 *
 *
 * To compile:
 * You need to put a copy of the nifti1.h header file in this directory.
 * It can be obtained from the NIFTI homepage  http://nifti.nimh.nih.gov/
 * or from the niftilib SourceForge site http://niftilib.sourceforge.net/
 *
 * cc -o nifti1_read_write nifti1_read_write.c
 *
 *
 * To run:
 * nifti1_read_write -w abc.nii abc.nii
 * nifti1_read_write -r abc.nii abc.nii
 *
 *
 * The read method is hardcoded to read float32 data.  To change
 * to your datatype, just change the line:
 * typedef float MY_DATATYPE;
 *
 * The write method is hardcoded to write float32 data.  To change
 * to your datatype, change the line:
 * typedef float MY_DATATYPE;
 * and change the lines:
 * hdr.datatype = NIFTI_TYPE_FLOAT32;
 * hdr.bitpix = 32;
 *
 *
 * Written by Kate Fissell, University of Pittsburgh, May 2005.
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <dirent.h>
#include <libgen.h>
#include <unistd.h>
#include <math.h>
#include "avs_io.h"
#include "nifti1.h"

#define MAX_PATH_LEN 1024
#define MY_DATATYPE float



#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

avs_header *nifti1_hdr_to_avs_hdr(nifti_1_header *hdr,
				  MY_DATATYPE minVal,
				  MY_DATATYPE maxVal){
  avs_header *newhdr = calloc(sizeof(avs_header),1);

  if (newhdr == NULL)
    return NULL;

  newhdr->ndim = hdr->dim[0];
  newhdr->dim1 = hdr->dim[1];
  newhdr->dim2 = hdr->dim[2];
  newhdr->dim3 = hdr->dim[3];
  newhdr->min_x = 0;
  newhdr->min_y = 0;
  newhdr->min_z = 0;
  newhdr->max_x = hdr->dim[1]*hdr->pixdim[1];
  newhdr->max_y = hdr->dim[2]*hdr->pixdim[2];
  newhdr->max_z = hdr->dim[3]*hdr->pixdim[3];
  newhdr->filetype = 0;
  newhdr->skip = 0;
  newhdr->nspace = 3;
  newhdr->veclen = 1;
  newhdr->dataname[0] = '\0';
 
  if ((minVal >= 0) && (maxVal <= 255))
    newhdr->datatype = 1;
  else
    newhdr->datatype = 3;

  return newhdr;
}

/**********************************************************************
 *
 * read_nifti_file
 *
 **********************************************************************/
double read_nifti_files(char *hdr_file, char *hdr_file2,char *data_file,
			       char *data_file2,MY_DATATYPE scaling){
nifti_1_header hdr;
nifti_1_header hdr2;
FILE *fp;
FILE *fp2;
int ret,ret2,i,j;
double tp,tn,fpp,fn;
double precision, recall,fscore,accuracy;
double total,total2;
MY_DATATYPE min, max;
MY_DATATYPE *data=NULL;
MY_DATATYPE *data2=NULL;
MY_DATATYPE min2,max2;
 unsigned char *bytebuffer;
 short *shortbuffer;


/********** open and read header */
fp = fopen(hdr_file,"r");
if (fp == NULL) {
        fprintf(stderr, "\nError opening header file %s\n",hdr_file);
        exit(1);
}
ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
if (ret != 1) {
        fprintf(stderr, "\nError reading header file %s\n",hdr_file);
        exit(1);
}
fclose(fp);


/********** print a little header information */
fprintf(stderr, "\n%s header information:",hdr_file);
fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix);
fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);
fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr.vox_offset));
fprintf(stderr, "\n");


/********** open the datafile, jump to data offset */
fp = fopen(data_file,"r");
if (fp == NULL) {
        fprintf(stderr, "\nError opening data file %s\n",data_file);
        exit(1);
}

ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);
if (ret != 0) {
        fprintf(stderr, "\nError doing fseek() to %ld in data file %s\n",(long)(hdr.vox_offset), data_file);
        exit(1);
}


/********** allocate buffer and read first 3D volume from data file */
data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
if (data == NULL) {
        fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
        exit(1);
}
ret = fread(data, sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
        fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",data_file,ret);
        exit(1);
}
fclose(fp);


/********** scale the data buffer  */
if (hdr.scl_slope != 0) {
        for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
                data[i] = (data[i] * hdr.scl_slope) + hdr.scl_inter;
}


/********** print mean of data */
 total = data[0];
 min = data[0];
 max = data[0]; 
 for (i=1; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++){
   total += data[i];
   if (data[i]<min) min = data[i];
   if (data[i]>max) max = data[i];
 }
 total /= (hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
 fprintf(stderr, "\nMean of ground truth volume in %s is %.3f\n",data_file,total);
 fprintf(stderr, "\nMax of ground truth volume in %s is %.3f\n",data_file,max);
 fprintf(stderr, "\nMin of  ground truth volume in %s is %.3f\n",data_file,min);

/**************************************************************************************************/

 /********** open and read header2 */
fp2 = fopen(hdr_file2,"r");
if (fp2 == NULL) {
        fprintf(stderr, "\nError opening header file %s\n",hdr_file2);
        exit(1);
}
ret2 = fread(&hdr2, MIN_HEADER_SIZE, 1, fp);
if (ret2!= 1) {
        fprintf(stderr, "\nError reading header file %s\n",hdr_file2);
        exit(1);
}
fclose(fp2);


/********** print a little header information */
fprintf(stderr, "\n%s header information:",hdr_file2);
fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr2.dim[1],hdr2.dim[2],hdr2.dim[3],hdr2.dim[4]);
fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr2.datatype,hdr2.bitpix);
fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr2.scl_slope,hdr2.scl_inter);
fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr2.vox_offset));
fprintf(stderr, "\n");


/********** open the datafile, jump to data offset */
fp2 = fopen(data_file2,"r");
if (fp2 == NULL) {
        fprintf(stderr, "\nError opening data file %s\n",data_file2);
        exit(1);
}

ret2 = fseek(fp2, (long)(hdr2.vox_offset), SEEK_SET);
if (ret2 != 0) {
        fprintf(stderr, "\nError doing fseek() to %ld in data file %s\n",(long)(hdr2.vox_offset), data_file2);
        exit(1);
}


/********** allocate buffer and read first 3D volume from data file */
data2 = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr2.dim[1]*hdr2.dim[2]*hdr2.dim[3]);
if (data2 == NULL) {
        fprintf(stderr, "\nError allocating data buffer for %s\n",data_file2);
        exit(1);
}
ret2 = fread(data2, sizeof(MY_DATATYPE), hdr2.dim[1]*hdr2.dim[2]*hdr2.dim[3], fp2);
if (ret2 != hdr2.dim[1]*hdr2.dim[2]*hdr2.dim[3]) {
        fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",data_file2,ret2);
        exit(1);
}
fclose(fp2);


/********** scale the data buffer  */
if (hdr2.scl_slope != 0) {
        for (i=0; i<hdr2.dim[1]*hdr2.dim[2]*hdr2.dim[3]; i++)
                data2[i] = (data2[i] * hdr2.scl_slope) + hdr2.scl_inter;
}


/********** print mean of data */
 total2 = data2[0];
 precision = data2[0];
 recall =data2[0];
 fscore = data2[0];
 accuracy = data2[0];
 tp=data2[0];
 tn=data2[0];
 fn=data2[0];
 fpp=data2[0];
 min2 = data2[0];
 max2 = data2[0]; 
 for (i=1; i<hdr2.dim[1]*hdr2.dim[2]*hdr2.dim[3]; i++){
   total2 += data[i];
   if (data2[i]<min2) min2 = data2[i];
   if (data2[i]>max2) max2 = data2[i];
 }
 total2 /= (hdr2.dim[1]*hdr2.dim[2]*hdr2.dim[3]);
 fprintf(stderr, "\nMean of classified volume in %s is %.3f\n",data_file2,total2);
 fprintf(stderr, "\nMax of classified volume in %s is %.3f\n",data_file2,max2);
 fprintf(stderr, "\nMin of classified volume in %s is %.3f\n",data_file2,min2);
 
for (i=1; i<hdr2.dim[1]*hdr2.dim[2]*hdr2.dim[3]; i++){
        if(data[i]==0  && data2[i]==0)tn+=1;
        else if (data[i]!=0&& data2[i]!=0) tp+=1;
        else if(data[i]!=0 && data2[i]==0)fn+=1;
        else if(data[i]==0 && data2[i]!=0)fpp+=1;}
precision = tp/(tp+fpp);
recall = tp/(tp+fn);
fscore = (2*precision*recall)/(precision+recall);
accuracy = (tp+tn)/(tp+tn+fn+fpp);
fprintf(stderr,"\n tp is %f\n",tp);
fprintf(stderr,"\n tn is %f\n",tn);
fprintf(stderr,"\n fn is %f\n",fn);
fprintf(stderr,"\nfp is %f\n",fpp);
fprintf(stderr,"\nprecision is %f\n",precision);
fprintf(stderr,"\nrecall is %f\n",recall);
fprintf(stderr,"\nfscore is %f\n",fscore);
fprintf(stderr,"\naccuracy is %f\n",accuracy);
return fscore;
}





int main(int argc, char *argv[]) {
    char *dir1_path = "/home/pooja/segtest"; // Change to the desired directory path
    char *dir2_path = "/home/pooja/outfiles90ts"; // Change to the desired directory path
    char path1[MAX_PATH_LEN];
    char path2[MAX_PATH_LEN];
    MY_DATATYPE scaling = 1.0;
    double num_pairs = 0.000000;
    double fscore = 0.000000;
    double fscore_avg=0.000000;
    double fscore_pair=0.000000;

    DIR *dir1 = opendir(dir1_path);
    if (dir1 == NULL) {
        fprintf(stderr, "Unable to open directory: %s\n", dir1_path);
        exit(1);
    }

    struct dirent *entry1;
    while ((entry1 = readdir(dir1)) != NULL) {
        if (entry1->d_type == DT_REG) { // regular file
            char *file_name = entry1->d_name; // get filename
            snprintf(path1, MAX_PATH_LEN, "%s/%s", dir1_path, file_name);
            snprintf(path2, MAX_PATH_LEN, "%s/%s", dir2_path, file_name);

            if(access(path2, F_OK) != -1) { // check if file exists in second directory
                fscore_pair = read_nifti_files(path1, path2, path1, path2, scaling);
                if (isnan(fscore_pair)) {
                fscore_pair = 0.0;}
                fscore += fscore_pair;
                num_pairs++;
            }
        }
    }
    closedir(dir1);

    if (num_pairs > 0) {
        fscore_avg = fscore / num_pairs;
        printf("Average fscore: %f\n", fscore_avg);
    } else {
        printf("No matching files found in both directories.\n");
    }

    exit(0);
}