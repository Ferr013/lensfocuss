#include <stdio.h>
#include "qfits.h"

int test_qfits() {
    printf("---------------------------\n");
    printf("QFITS library test program.\n");

    char * fname = "/Users/giofer/Documents/Programming/lensing_C/data/fits/image.fits";
    printf("Reading FITS file: %s\n", fname);

    int n_ext = qfits_query_n_ext(fname);
    if(n_ext == -1){
        printf("Failed to read header.\n");
        exit(1);
    }
    printf("Number of extensions: %d\n", n_ext);

    qfits_header* hdr = qfits_header_readext(fname, 0);
    if (!hdr) {
        printf("Failed to read header.\n");
        exit(1);
    }

    int seg_start, seg_size;
    int hdr_info = qfits_get_hdrinfo(fname, 0, &seg_start, &seg_size);
    if (hdr_info == -1) {
        printf("Failed to get header info.\n");
        exit(1);
    }

    qfitsloader *ql = malloc(sizeof(qfitsloader));
    if (!ql) {
        printf("Failed to allocate memory for QFITS loader.\n");
        exit(1);
    }
    ql->filename = fname;
    ql->xtnum = 0;
    ql->pnum = 0;
    if (qfitsloader_init(ql)!=0) {
        printf("cannot read info about %s\n", fname);
        return -1 ;
    }

    printf( "file         : %s\n"
            "xtnum        : %d\n"
            "pnum         : %d\n"
            "# xtensions  : %d\n"
            "size X       : %d\n"
            "size Y       : %d\n"
            "planes       : %d\n"
            "bitpix       : %d\n"
            "datastart    : %d\n"
            "datasize     : %d\n"
            "bscale       : %g\n"
            "bzero        : %g\n",
            ql->filename,
            ql->xtnum,
            ql->pnum,
            ql->exts,
            ql->lx,
            ql->ly,
            ql->np,
            ql->bitpix,
            ql->seg_start,
            ql->seg_size,
            ql->bscale,
            ql->bzero);

    if(qfits_loadpix(ql) != 0) {
        printf("Failed to load pixel data.\n");
        free(ql);
        return -1;
    }

    printf("---------------------------\n");
    for(int i = 0; i < ql->lx; i++) {
        for(int j = 0; j < ql->ly; j++) {
            if((i*j)%101 == 0) printf("Pixel (%d, %d): %f\n", i, j, ql->fbuf[i * ql->ly + j]);
        }
    }
    printf("---------------------------\n");
    free(ql);
    return 0;
}

// int main() {
//     test_qfits();
//     return 0;
// }