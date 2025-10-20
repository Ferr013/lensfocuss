/*******************************************************************************************
*
*   Copyright (c) 2025 Giovanni Ferrami (@Ferr013)
*
********************************************************************************************/
#include <stdlib.h>
#include "raylib.h"
#include "raymath.h"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#undef RAYGUI_IMPLEMENTATION // Avoid including raygui implementation again

#define GUI_WINDOW_FILE_DIALOG_IMPLEMENTATION
#include "gui_window_file_dialog.h"


#include "style_amber.h"
#include "qfits.h" //qfits version 6.20 by yjung
#include "fastell.h"

//------------------------------------------------------------------------------------
// Debugging
#define DEBUG 1
#define CURRENT_LOG_MODE LOG_INFO
#define READ_FITS_IMAGE 1
char img_fname[512] = { 0 };
const char default_img_fname[] = "data/fits/J0206_SCI.fits";
// strcpy(img_fname, default_img_fname); // At beginning of main()
//------------------------------------------------------------------------------------
// Images colors
bool DARKMODE = true;
const Color colors[] = {ORANGE, BLUE, RED, PURPLE, GOLD, LIME, DARKPURPLE, BEIGE};
const int numColors = sizeof(colors) / sizeof(colors[0]);
bool DRAW_EXT_AS_POINTS = true;
bool DRAW_EXT_PARITY    = false;
//------------------------------------------------------------------------------------
// Constants
#define DISTANCE_LENS_TO_OBSERVER 1.0         // Distance to the observer (arbitrary units)
#define DISTANCE_LENS_TO_SOURCE 2.0           // Distance from lens to the source (arbitrary units)
#define MASS_LENS 1.0                         // Mass of the lens (arbitrary units)
#define RAND_SRC_POS 2                        // Radius in arcsecs in which to spawn a source
#define RAND_LENS_POS 1                       // Radius in arcsecs in which to spawn a lens
#define R_SCHWARTZ 2953.339382                // Schwarzschild radius for 1 solar mass in meters
#define MPC2METER 3.085677581e22              // Mpc to meters
#define ARCSEC2RADIAN 0.00000484813681109536  // <---
#define BLACKBOARD CLITERAL(Color){ 34, 34, 34, 255} // Blackboard color
#define MAXIMUM_DYNARRAYSIZE 100000
#define PIUNICODE 0x03c0
//------------------------------------------------------------------------------------
//  TERMINAL COLORS
#define TERM_RED     "\033[31m"
#define TERM_GREEN   "\033[32m"
#define TERM_YELLOW  "\033[33m"
#define TERM_BLUE    "\033[34m"
#define TERM_MAGENTA "\033[35m"
#define TERM_CYAN    "\033[36m"
#define TERM_RESET   "\033[0m"
//------------------------------------------------------------------------------------
// Angular diameter distances in Mpc (hardcoded for now)
#define DIST_OL 1259.0251822341118 // zl = 0.5
#define DIST_OS 1588.5730743925365 // zs = 3
#define DIST_LS 1116.4386310547445 // zl = 0.5, zs = 3
//------------------------------------------------------------------------------------
const double LENS_MASS = 8e11;      // Mass of the main deflector
const int    N_EXTSRC  = 512;       // n point to discretize circular extended source

const int    N_IGRIDP       = 80 ;    // N of points in img plane grid [NxN]
const int    I_to_S_SPACING = 3;      // sampling every 2x2 pixel in img plane to
                                      // create the src plane
const int    SRC_PLANE_SIZE = 50;     // N of points in src plane grid [NxN]

double IMG_WIDTH_AS  = 8.0;   // Width and ...
double IMG_HEIGHT_AS = 8.0;   // Heigth of the image plane in arcsec

double       R_EXTSRC  = 0.06;      // radius of the extended source, in arcsec
double       angStep   = 2*PI/(double)N_EXTSRC;
double   SCALE_SPLANE  = 0.01;      // Scale of 1 pixel in the source plane in arcsecs
double   SCALE_IPLANE  = 0.02;      // Scale of 1 pixel in the source plane in arcsecs
double       GammaExt1 = 0;         // First component of the external shear
double       GammaExt2 = 0;         // Second component of the external shear

//------------------------------------------------------------------------------------
const int screenWidth = 1512;
const int screenHeight = 912;

const int sizePlanes = 600;
const int paddingUIx  = 30;
const int paddingUIy  = 10;

const int panelUIposx = 2 * sizePlanes + paddingUIx / 2;
const int panelUIposy = paddingUIy;
const int panelUIsizx = 270;
const int panelUIsizy = screenHeight - paddingUIy * 2;

const int widthExtraP  = screenWidth - panelUIsizx - paddingUIx - 5;
const int heigthExtraP = screenHeight - sizePlanes - paddingUIy*2;

const int ctrSrcPlaneX = (screenWidth - panelUIsizx - paddingUIx) / 4;
const int ctrSrcPlaneY = sizePlanes / 2;
const int ctrImgPlaneX = (screenWidth - panelUIsizx - paddingUIx) / 4 * 3;
const int ctrImgPlaneY = sizePlanes / 2;
const int ctrExtraPanX = (screenWidth - panelUIsizx - paddingUIx) / 2;
const int ctrExtraPanY = screenHeight - (screenHeight - sizePlanes - paddingUIy) / 2 - 10;

const int zeroSrcPlaneX = ctrSrcPlaneX - sizePlanes / 2;
const int zeroSrcPlaneY = ctrSrcPlaneY - sizePlanes / 2;
const int zeroImgPlaneX = ctrImgPlaneX - sizePlanes / 2;
const int zeroImgPlaneY = ctrImgPlaneY - sizePlanes / 2;
const int zeroExtraPanX = ctrExtraPanX - widthExtraP / 2;
const int zeroExtraPanY = ctrExtraPanY - heigthExtraP/ 2;

Rectangle srcPlaneBounds = (Rectangle){zeroSrcPlaneX, zeroSrcPlaneY, sizePlanes, sizePlanes};
Rectangle imgPlaneBounds = (Rectangle){zeroImgPlaneX, zeroImgPlaneY, sizePlanes, sizePlanes};
Rectangle otrPlaneBounds = (Rectangle){zeroExtraPanX, zeroExtraPanY, widthExtraP, heigthExtraP};

double gridLimDistance = 10;

//------------------------------------------------------------------------------------
enum MODE { // Toggle between behaviors
    LENS_TO_SRC_MODE,
    SRC_TO_LENS_MODE,
    NUM_MODES,
  } CURRENT_MODE;

// UI elements to switch mode
const char* Get_Mode_Text(enum MODE CURRENT_MODE) {
    static char* replies[] = {
        "Lens --> Source ",
        "Lens <-- Source ",
        "NUM_MODES"};
        return replies[CURRENT_MODE];
}

// **QFITS wrapper**
//------------------------------------------------------------------------------------
typedef struct {
    char * filename;
    int lx;
    int ly;
    float  *fBuffer;
    double *dBuffer;
    bool isDouble;
} QFitsInterface;

QFitsInterface* ReadFitsFile(char* fname){
    printf("Reading FITS file: %s\n", fname);

    qfits_header* hdr = qfits_header_readext(fname, 0);
    if (!hdr) {
        TraceLog(LOG_ERROR, "Failed to read header.\n");
        exit(1);
    }

    int seg_start, seg_size;
    int hdr_info = qfits_get_hdrinfo(fname, 0, &seg_start, &seg_size);
    if (hdr_info == -1) {
        TraceLog(LOG_ERROR, "Failed to get header info.\n");
        exit(1);
    }

    qfitsloader *ql = malloc(sizeof(qfitsloader));
    if (!ql) {
        TraceLog(LOG_ERROR, "Failed to allocate memory for QFITS loader.\n");
        exit(1);
    }
    ql->filename = fname;
    ql->xtnum = 0;
    ql->pnum = 0;
    if (qfitsloader_init(ql)!=0) {
        TraceLog(LOG_ERROR, "cannot read info about %s\n", fname);
        exit(1);
    }

    TraceLog(LOG_INFO,
            "file         : %s\n"
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
        TraceLog(LOG_ERROR, "Failed to load pixel data.\n");
        free(ql);
        exit(1);
    }

    QFitsInterface* fitsData = malloc(sizeof(QFitsInterface));
    fitsData->filename = ql->filename;
    fitsData->lx = ql->lx;
    fitsData->ly = ql->ly;

    ql->bitpix = -32; //REMOVE!

    if(ql->bitpix == -32) { // float
        fitsData->isDouble = false;
        fitsData->fBuffer = (float*)calloc(fitsData->lx * fitsData->ly, sizeof(float));
        if (!fitsData->fBuffer) {
            TraceLog(LOG_ERROR, "Failed to allocate memory for FITS buffer.\n");
            free(fitsData);
            free(ql);
            return NULL;
        }
        TraceLog(LOG_DEBUG, "Allocated memory for FITS float buffer.\n");
        memcpy(fitsData->fBuffer, ql->fbuf, fitsData->lx * fitsData->ly * sizeof(float));

    } else if(ql->bitpix == -64) { // double
        fitsData->isDouble = true;
        fitsData->dBuffer = (double*)calloc(fitsData->lx * fitsData->ly, sizeof(double));
        if (!fitsData->dBuffer) {
            TraceLog(LOG_ERROR, "Failed to allocate memory for FITS buffer.\n");
            free(fitsData);
            free(ql);
            return NULL;
        }
        TraceLog(LOG_DEBUG, "Allocated memory for FITS double buffer.\n");
        memcpy(fitsData->dBuffer, ql->dbuf, fitsData->lx * fitsData->ly * sizeof(double));
    }
    TraceLog(LOG_DEBUG, "Successfully allocated memory for FITS buffer.\n");
    free(ql);
    return fitsData;
}

void ClearFitsData(QFitsInterface* fitsData) {
    if (fitsData) {
        free(fitsData->fBuffer);
        free(fitsData->dBuffer);
        free(fitsData);
    }
}

// **Math**
//------------------------------------------------------------------------------------
typedef double (*func_ptr)(double);

typedef struct {
    double x;
    double y;
} Vector2Double;

int sign(double x) {
    return (x > 0) - (x < 0);
}

double SIGN(double x, double y) {
    return (y >= 0.0) ? fabs(x) : -fabs(x);
}

double RootFindBrent(func_ptr func, double x1, double x2, double tol) { //From Numerical Recipes - Sect 9.3
    const int ITMAX = 100;
    const double EPS = 1e-5;
    double a=x1,b=x2,c=x2,d,e,fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        return -9999.;
    fc=fb;
    for (int iter=0;iter<ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) {
            return b;
        }
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=fabs(p);
            double min1=3.0*xm*q-fabs(tol1*q);
            double min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
            b += d;
        else {
            b += SIGN(tol1,xm);
            fb=func(b);
        }
    }
    return -9999.;
}
// **Utils**
//------------------------------------------------------------------------------------
double RandOnUnitBall() { // Random value between -1 and 1
    return ((double)rand() / RAND_MAX) * 2.0 - 1.0;
}

void physics_util_translate_rotate_fwd(double x, double y, double x0, double y0,
		double sinth, double costh, double *x_out, double *y_out)
{
	double sx, sy;

	// translate
	sx = x - x0;
	sy = y - y0;

	// rotate
	*x_out = sx*costh + sy*sinth;
	*y_out = -sx*sinth + sy*costh;
}


void physics_util_rotate_bwd(double x, double y,
		double sinth, double costh, double *x_out, double *y_out)
{
	// inverse rotation
	*x_out = x*costh - y*sinth;
	*y_out = x*sinth + y*costh;
}

// **Draw pixels**
//------------------------------------------------------------------------------------
Color getColorFromValue(double value){
    float r, g, b, a;
    float r1,g1,b1;
    float r2,g2,b2;
    // Clamp value between 0 and 1
    if (value > 1.0f) value = 1.0f;
    if (value < -1.0f) value = -1.0f;

    // A colormap that goes from green to purple to orange
    //                          [(-1.0      0.0       1.0]
    if (value < 0.0f) {
        value = -value; // Use absolute value for color mapping
        // Define purple (at value = 0.0)
        r1 = 50.0f;
        g1 = 0.0f;
        b1 = 75.0f;

        // Define orange (at value = (-)1.0)
        r2 = 0.0f;
        g2 = 165.0f;
        b2 = 75.0f;
    } else {
        // Define purple (at value = 0.0)
        r1 = 50.0f;
        g1 = 0.0f;
        b1 = 75.0f;

        // Define orange (at value = 1.0)
        r2 = 255.0f;
        g2 = 165.0f;
        b2 = 0.0f;
    }
    // Linear interpolation between purple and orange

    if (value == 0.0f) {
        r = 100; g = 100; b = 100; a = 255.0;
    } else {
        r = r1 + (r2 - r1) * value;
        g = g1 + (g2 - g1) * value;
        b = b1 + (b2 - b1) * value;
        a = 255.0;
    }

    Color pixelColor = {
        (unsigned char)r,
        (unsigned char)g,
        (unsigned char)b,
        (unsigned char)a  // Full opacity
    };
    return pixelColor;
}
void UpdateTextureLensfocus(Texture2D* texture, double* img, int width, int height){
    Image canvas = GenImageColor(texture->width, texture->height, BLANK);

    for (int n = 0; n < texture->width*texture->height; n++) {
        Color pixelColor = getColorFromValue(img[n]);
        ImageDrawPixel(&canvas, n % width, n / width, pixelColor);
    }
    *texture = LoadTextureFromImage(canvas);
    UnloadImage(canvas);
}
// **Coordinates**
//------------------------------------------------------------------------------------
typedef struct{
    int x;
    int y;
} PixelCoord;

typedef struct{
    float RA;
    float DEC;
} SkyCoord;

typedef struct {
    SkyCoord images[N_EXTSRC];
} ExtSource; //TODO: Make this a DynamicArraySkyCoord

// Dynamic array structure
typedef struct {
    SkyCoord* coords;
    size_t size;
    size_t capacity;
} DynamicArraySkyCoord;

void InitDynArraySkyCoord(DynamicArraySkyCoord* array, size_t initialCapacity) {
    array->coords = malloc(initialCapacity * sizeof(SkyCoord));
    if (array->coords == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    array->size = 0;
    array->capacity = initialCapacity;
}

void AppendItemDynArraySkyCoord(DynamicArraySkyCoord* array, SkyCoord item) {
    if (array->size >= MAXIMUM_DYNARRAYSIZE){
        TraceLog(LOG_ERROR, "Maximum size reached\n");
        exit(1);
    }
    if (array->size >= array->capacity) {
        size_t newCapacity = array->capacity * 2;
        SkyCoord* newCoords = realloc(array->coords, newCapacity * sizeof(SkyCoord));
        if (newCoords == NULL) {
            TraceLog(LOG_ERROR, "Memory reallocation failed\n");
            exit(1);
        }
        array->coords = newCoords;
        array->capacity = newCapacity;
        TraceLog(LOG_DEBUG, "Array resized to capacity: %zu\n", array->capacity);
    }
    array->coords[array->size] = item;
    array->size++;
}

void CleanDynArraySkyCoord(DynamicArraySkyCoord* array, size_t initialCapacity) {
    free(array->coords);
    InitDynArraySkyCoord(array, initialCapacity);
}

void FreeDynArraySkyCoord(DynamicArraySkyCoord* array) {
    free(array->coords);
    array->coords = NULL;
    array->size = 0;
    array->capacity = 0;
}

typedef struct {
    DynamicArraySkyCoord* imgSCoordDArray;
    size_t size;
    size_t capacity;
} DynamicArrayImages;

void InitDAImages(DynamicArrayImages* array, size_t initialCapacity) {
    array->imgSCoordDArray = malloc(initialCapacity * sizeof(DynamicArraySkyCoord));
    if (array->imgSCoordDArray == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    array->size = 0;
    array->capacity = initialCapacity;
}

void AppendItemDAImages(DynamicArrayImages* array, DynamicArraySkyCoord item) {
    if (array->size >= MAXIMUM_DYNARRAYSIZE){
        TraceLog(LOG_ERROR, "Maximum size reached\n");
        exit(1);
    }
    if (array->size >= array->capacity) {
        size_t newCapacity = array->capacity * 2;
        DynamicArraySkyCoord* newImgArray = realloc(array->imgSCoordDArray, newCapacity * sizeof(DynamicArraySkyCoord));
        if (newImgArray == NULL) {
            TraceLog(LOG_ERROR, "Memory reallocation failed\n");
            exit(1);
        }
        array->imgSCoordDArray = newImgArray;
        array->capacity = newCapacity;
        TraceLog(LOG_DEBUG, "Array resized to capacity: %zu\n", array->capacity);
    }
    array->imgSCoordDArray[array->size] = item;
    array->size++;
}

void CleanDAImages(DynamicArrayImages* array, size_t initialCapacity) {
    free(array->imgSCoordDArray);
    InitDAImages(array, initialCapacity);
}

void FreeDAImages(DynamicArrayImages* array) {
    free(array->imgSCoordDArray);
    array->imgSCoordDArray = NULL;
    array->size = 0;
    array->capacity = 0;
}

PixelCoord MapSkyToPixel(SkyCoord skyCoord, int screenWidth, int screenHeight, bool SRCPLANE) {
    PixelCoord pixelCoord;
    if (SRCPLANE == 1) {
        pixelCoord.x = ctrSrcPlaneX + (int)(skyCoord.RA  / SCALE_SPLANE);
        pixelCoord.y = ctrSrcPlaneY + (int)(skyCoord.DEC / SCALE_SPLANE);
    } else {
        pixelCoord.x = ctrImgPlaneX + (int)(skyCoord.RA  / SCALE_IPLANE);
        pixelCoord.y = ctrImgPlaneY + (int)(skyCoord.DEC / SCALE_IPLANE);
    }
    return pixelCoord;
}

SkyCoord MapPixeltoSky(PixelCoord pixelCoord, int screenWidth, int screenHeight, bool SRCPLANE) {
    SkyCoord skyCoord;
    if (SRCPLANE == 1) {
        skyCoord.RA = (pixelCoord.x - ctrSrcPlaneX)*SCALE_SPLANE;
        skyCoord.DEC = (pixelCoord.y - ctrSrcPlaneY)*SCALE_SPLANE;
    } else {
        skyCoord.RA = (pixelCoord.x - ctrImgPlaneX)*SCALE_IPLANE;
        skyCoord.DEC = (pixelCoord.y - ctrImgPlaneY)*SCALE_IPLANE;
    }
    return skyCoord;
}

SkyCoord MapSkyPolarToSkyCartesian (double r, double theta) {
    SkyCoord skyCoordCart;
    skyCoordCart.RA  = r * sin(theta);
    skyCoordCart.DEC = r * cos(theta);
    return skyCoordCart;
}

// ** Draw interface
//------------------------------------------------------------------------------------

int CheckPixelInTargetWindow(PixelCoord p, int TargetWindow) {
    if (TargetWindow == 0) { // Image plane
        if (p.x > imgPlaneBounds.x &&
            p.x < imgPlaneBounds.x + imgPlaneBounds.width &&
            p.y > imgPlaneBounds.y &&
            p.y < imgPlaneBounds.y + imgPlaneBounds.height
            #if 1 // Keeping the window label area clear
            && !(p.x > imgPlaneBounds.x + 10  &&
            p.x < imgPlaneBounds.x + 10 + 130 &&
            p.y > imgPlaneBounds.y + 15 &&
            p.y < imgPlaneBounds.y + 15 + 25)
#endif
            )
            {
            return 1;
        } else { return 0; }
    }
    if (TargetWindow == 1) { // Source plane
        if (
            p.x > srcPlaneBounds.x &&
            p.x < srcPlaneBounds.x + srcPlaneBounds.width &&
            p.y > srcPlaneBounds.y &&
            p.y < srcPlaneBounds.y + srcPlaneBounds.height
#if 1 // Keeping the window label area clear
            && !(p.x > srcPlaneBounds.x + 5  &&
            p.x < srcPlaneBounds.x + 5 + 145 &&
            p.y > srcPlaneBounds.y + 10 &&
            p.y < srcPlaneBounds.y + 10 + 30)
#endif
            )
        {
            return 1;
        } else { return 0; }
    }
    if (TargetWindow == 2) { // Other stuff window
        if (p.x > otrPlaneBounds.x &&
            p.x < otrPlaneBounds.x + otrPlaneBounds.width &&
            p.y > otrPlaneBounds.y &&
            p.y < otrPlaneBounds.y + otrPlaneBounds.height
#if 1 // Keeping the window label area clear
            && !(p.x > otrPlaneBounds.x + 5  &&
            p.x < otrPlaneBounds.x + 5 + 145 &&
            p.y > otrPlaneBounds.y + 10 &&
            p.y < otrPlaneBounds.y + 10 + 30)
#endif
            )
            {
            return 1;
        } else { return 0; }
    }
    return -1; // Invalid target window
}

void SetupSrcImgPlanes() {
    ClearBackground(DARKMODE ? BLACKBOARD : RAYWHITE);
    DrawRectangle(zeroSrcPlaneX, zeroSrcPlaneY, sizePlanes, sizePlanes, DARKMODE ? BLACKBOARD : RAYWHITE);
    DrawRectangleLinesEx(srcPlaneBounds, 5, RED);
    DrawText("Source plane", 15, 15, 20, RED);
    DrawRectangleLinesEx(imgPlaneBounds, 5, DARKMODE ? SKYBLUE : BLUE);
    DrawText("Image plane", zeroImgPlaneX + 15, 15, 20, DARKMODE ? SKYBLUE : BLUE);
    DrawRectangleLinesEx(otrPlaneBounds, 5, DARKMODE ? RAYWHITE : BLACKBOARD);
    // DrawText("Other stuff", zeroExtraPanX + 15, zeroExtraPanY + 15, 20, DARKMODE ? RAYWHITE : BLACKBOARD);
}

void DrawScaleBar() {
    SkyCoord oneArcsec = {1.0, 0.0};
    int widthSRC = (int)(oneArcsec.RA  / SCALE_SPLANE);
    int widthIMG = (int)(oneArcsec.RA  / SCALE_IPLANE);

    int height = 8;
    Color color = DARKMODE ? RAYWHITE : BLACK;

    DrawText("1 \"", zeroSrcPlaneX + sizePlanes - 120, zeroSrcPlaneY + sizePlanes - 30 - height - 9, 20, color);
    DrawRectangle(zeroSrcPlaneX + sizePlanes - 110 - widthSRC/2, zeroSrcPlaneY + sizePlanes - 30, widthSRC, height, color);
    DrawRectangleLinesEx((Rectangle){zeroSrcPlaneX + sizePlanes - 110 - widthSRC/2, zeroSrcPlaneY + sizePlanes - 30, widthSRC, height}, 2, DARKGRAY);

    DrawText("1 \"", zeroImgPlaneX + sizePlanes - 120, zeroImgPlaneY + sizePlanes - 30 - height - 9, 20, color);
    DrawRectangle(zeroImgPlaneX + sizePlanes - 110 - widthIMG/2, zeroImgPlaneY + sizePlanes - 30, widthIMG, height, color);
    DrawRectangleLinesEx((Rectangle){zeroImgPlaneX + sizePlanes - 110 - widthIMG/2, zeroImgPlaneY + sizePlanes - 30, widthIMG, height}, 2, DARKGRAY);
}

//------------------------------------------------------------------------------------

typedef struct {
    double D_OL;
    double D_OS;
    double D_LS;
    double lens_dist_ratio;
} AngCosmoDistances;

// **Source and Image Planes**
// This are adapted from G. Vernardos 'veryknotty' code.
//------------------------------------------------------------------------------------
typedef struct  {
    int Ni;                        //pixels in x direction
    int Nj;                        //pixels in y direction
    int Nm;                        //total pixels in the image data
    double width, height;          //in arcsec
    double xmin, xmax, ymin, ymax;
    double* img;                   //values of the pixels
    double* x;                     //pixel x-coordinates in arcsec
    double* y;                     //pixel y-coordinates in arcsec
    double* defl_x;                //deflected x coordinates
    double* defl_y;                //deflected y coordinates

    QFitsInterface* fitsImageData;
    QFitsInterface* fitsPSFData;

    Texture2D* textureImageData;
 } ImagePlane;

void LoadFitsFileInImagePlane(ImagePlane* imgPlane, char* img_fname, double img_width_as, double img_height_as){
    imgPlane->fitsImageData = ReadFitsFile(img_fname);

    if (!imgPlane->fitsImageData) {
        TraceLog(LOG_ERROR, "Failed to create FITS interface.\n");
        imgPlane->Ni     = N_IGRIDP;                    //pixels in x direction
        imgPlane->Nj     = N_IGRIDP;                    //pixels in y direction
        imgPlane->Nm     = imgPlane->Ni * imgPlane->Nj; //total pixels in the image data
        imgPlane->width  = img_width_as;                //in arcsec
        imgPlane->height = img_height_as;               //in arcsec
        imgPlane->img     = (double*) calloc(imgPlane->Nm,sizeof(double));
    } else {
        imgPlane->Ni     = imgPlane->fitsImageData->lx; //pixels in x direction
        imgPlane->Nj     = imgPlane->fitsImageData->ly; //pixels in y direction
        imgPlane->Nm     = imgPlane->Ni * imgPlane->Nj; //total pixels in the image data
        // TODO: read the size in arcsecs from the fits file!
        imgPlane->width  = img_width_as;                //in arcsec
        imgPlane->height = img_height_as;               //in arcsec
        imgPlane->img     = (double*) calloc(imgPlane->Nm,sizeof(double));

        for(int i=0;i<imgPlane->Ni;i++){
            for(int j=0;j<imgPlane->Nj;j++){
                //convert FITS standard (bottom to top) to the one used in this code (top to bottom)
                if (imgPlane->fitsImageData->isDouble) {
                    imgPlane->img[i*imgPlane->Nj+j] = imgPlane->fitsImageData->dBuffer[(imgPlane->fitsImageData->ly-i-1)*imgPlane->fitsImageData->lx+j];
                } else {
                    imgPlane->img[i*imgPlane->Nj+j] = (double)imgPlane->fitsImageData->fBuffer[(imgPlane->fitsImageData->ly-i-1)*imgPlane->fitsImageData->lx+j];
                }
            }
        }
    }
    imgPlane->textureImageData = (Texture2D*) malloc(sizeof(Texture2D));
    TraceLog(LOG_DEBUG, "Updating texture from image data");
    imgPlane->textureImageData->width  = imgPlane->Ni;
    imgPlane->textureImageData->height = imgPlane->Nj;
    UpdateTextureLensfocus(imgPlane->textureImageData, imgPlane->img, imgPlane->Ni, imgPlane->Nj);
    TraceLog(LOG_DEBUG, "Finished updating texture from image data");
}

void InitImagePlane(ImagePlane* imgPlane, char* img_fname, double img_width_as, double img_height_as){
    LoadFitsFileInImagePlane(imgPlane, img_fname, img_width_as, img_height_as);
    imgPlane->x       = (double*) calloc(imgPlane->Nm,sizeof(double));
    imgPlane->y       = (double*) calloc(imgPlane->Nm,sizeof(double));
    imgPlane->defl_x  = (double*) calloc(imgPlane->Nm,sizeof(double));
    imgPlane->defl_y  = (double*) calloc(imgPlane->Nm,sizeof(double));

    int i0    = floor(imgPlane->Ni/2);
    int j0    = floor(imgPlane->Nj/2);
    double di = imgPlane->height/imgPlane->Ni;
    double dj = imgPlane->width/imgPlane->Nj;

    // xmax,xmin are the x coordinates of the leftmost and rightmost pixels.
    // REMEMBER: xmax-xmin != width (Nj*dj).
    // Similarly for y.
    imgPlane->xmin = -j0*dj;
    imgPlane->xmax = (imgPlane->Nj-1-j0)*dj;
    imgPlane->ymin = -i0*di;
    imgPlane->ymax = (imgPlane->Ni-1-i0)*di;

    for(int i=0;i<imgPlane->Ni;i++){
        for(int j=0;j<imgPlane->Nj;j++){
            imgPlane->x[i*imgPlane->Nj+j] =  (j-j0)*dj;
            imgPlane->y[i*imgPlane->Nj+j] = -(i-i0)*di;//reflect y-axis
        }
    }
}

void UpdateImagePlane(ImagePlane* imgPlane, char* img_fname, double img_width_as, double img_height_as){
    free(imgPlane->img);
    free(imgPlane->x);
    free(imgPlane->y);
    free(imgPlane->defl_x);
    free(imgPlane->defl_y);
    ClearFitsData(imgPlane->fitsImageData);
    UnloadTexture(*imgPlane->textureImageData);
    free(imgPlane->textureImageData);
    InitImagePlane(imgPlane, img_fname, img_width_as, img_height_as);
}

// **Source and Image positions**
//------------------------------------------------------------------------------------
SkyCoord RandomSrcPosition() {
    SkyCoord pcoord = {RandOnUnitBall()*RAND_SRC_POS, RandOnUnitBall()*RAND_SRC_POS};
    return pcoord;
}

SkyCoord RandomLensPosition() {
    SkyCoord pcoord = {RandOnUnitBall()*RAND_LENS_POS, RandOnUnitBall()*RAND_LENS_POS};
    return pcoord;
}

// **Lens Manager**
//------------------------------------------------------------------------------------

typedef enum {
    LTYPE_POINT,
    LTYPE_SIE,
    LTYPE_SPEMD,
    LTYPE_PERT
} LensModelType;

LensModelType SELECTED_LENS_MODEL = LTYPE_SPEMD;
//------------------------------------------------------------------------------------
// NOTE: Dropdown menu require GuiLock(), so we set them here
bool dropdownLensModelEditMode = false;
int  dropdownLensModel = 2;
const char *dropdownLensModelItems = "Point;SIE;SPEMD";

bool dropdownSourceModelEditMode = false;
int  dropdownSourceModel = 1;
const char *dropdownSourceModelItems = "Point;Extended";

typedef struct
{
    double M;
    SkyCoord coord;       // Lens position [arcsec]
    LensModelType modelType;

    double b;              // Could be Einstein radius [arcsec]
    double qh;             // proj (?) slope : \gamma = 2*qh+1
    double rc;             // core radius [arcsec]
    double f, PA;          // Axis ratio and Position angle [radians]
    double cos_PA, sin_PA; // precomputed for efficiency
} Lens;

typedef struct {
    Lens* lenses;
    int count;
    int capacity;
} LensManager;

void InitLensManager(LensManager* manager){
    manager->capacity = 10;
    manager->count    = 0;
    manager->lenses   = malloc(sizeof(Lens) * manager->capacity);
}

int AddLens(LensManager* manager, Lens* lens) {
    if (manager->count >= manager->capacity) {
        manager->capacity *= 2;
        manager->lenses = realloc(manager->lenses, sizeof(Lens) * manager->capacity);
    }
    manager->lenses[manager->count] = *lens;
    return manager->count++;
}

void RemoveLens(LensManager* manager, int index) {
    if (index < 0 || index >= manager->count) return;
    manager->lenses[index] = manager->lenses[manager->count - 1];
    manager->count--;
}

void FreeLensManager(LensManager* manager) {
    free(manager->lenses);
    manager->lenses = NULL;
    manager->count = manager->capacity = 0;
    free(manager);
}

Lens* GetRandomLens(){
    Lens* lens = malloc(sizeof(Lens));
    lens->M = (1.25 + 0.75 * RandOnUnitBall()) * 1e12;
    lens->coord = (SkyCoord){0., 0.};
    lens->modelType = SELECTED_LENS_MODEL;
    lens->b = (1.0 + 0.5 * RandOnUnitBall());
    lens->rc = 1e-6;
    lens->qh = (0.5 + 0.2 * RandOnUnitBall());
    lens->f  = 0.2+(RandOnUnitBall()+1.0)/2.0*0.79;
    lens->PA = PI*(0.5 + RandOnUnitBall()/2.0);
    lens->cos_PA = cos(lens->PA);
    lens->sin_PA = sin(lens->PA);
    lens->b  = 1.25+(RandOnUnitBall()+1.0)/1.5;
    return lens;
}

void ReplaceLastLensWithRandomLens(LensManager* manager) {
    if (manager->count>0){
        RemoveLens(manager, manager->count - 1);
    }
    Lens* lens = GetRandomLens();
    int lens_id = AddLens(manager, lens);
    free(lens);
}

// **Source Manager**
//------------------------------------------------------------------------------------
typedef enum {
    STYPE_POINT,
    STYPE_EXT
} SourceModelType;

typedef struct
{
  SkyCoord coord;
  SourceModelType modelType;
  double radius;
} Source;

typedef struct {
    Source* sources;
    int count;
    int capacity;
} SourceManager;

void InitSourceManager(SourceManager* manager){
    manager->capacity = 10;
    manager->count    = 0;
    manager->sources  = malloc(sizeof(Source) * manager->capacity);
}

int AddSource(SourceManager* manager, Source* source) {
    if (manager->count >= manager->capacity) {
        manager->capacity *= 2;
        manager->sources = realloc(manager->sources, sizeof(Source) * manager->capacity);
    }
    manager->sources[manager->count] = *source;
    return manager->count++;
}

void RemoveSource(SourceManager* manager, int index) {
    if (index < 0 || index >= manager->count) return;
    manager->sources[index] = manager->sources[manager->count - 1];
    manager->count--;
}

void FreeSourceManager(SourceManager* manager) {
    free(manager->sources);
    manager->sources = NULL;
    manager->count = manager->capacity = 0;
    free(manager);
}

Source* GetRandomSource(){
    Source* source = malloc(sizeof(Source));
    source->coord = RandomSrcPosition();
    source->modelType = STYPE_EXT;
    source->radius = R_EXTSRC;
    return source;
}

void ReplaceLastSourceWithRandomSource(SourceManager* manager) {
    if (manager->count>0){
        RemoveSource(manager, manager->count - 1);
    }
    Source* source = GetRandomSource();
    int source_id = AddSource(manager, source);
    free(source);
}


// **State**
//------------------------------------------------------------------------------------
// TODO: Refactor everything with arena allocators
typedef struct {
    LensManager* lensMan;
    SourceManager* sourceMan;

    SkyCoord obsImg[8]; // Up to 8 image positions for point sources
    bool obsImgIsVisible[8]; // Visibility flags for each image position

    DynamicArrayImages* images;
    DynamicArrayImages* caustics;
    DynamicArrayImages* criticalCurves;

    PixelCoord ctrSrcPlane; // I guess this could change, boh
    PixelCoord ctrImgPlane;

    int grid_size;
    DynamicArraySkyCoord* imgSGrid;
    DynamicArraySkyCoord* srcSGrid;

    ImagePlane* imgPlane;

    float gamma1;
    float gamma2;
} State;

State* StateInit () {
    State* state = (State*)malloc(sizeof(State));
    state->lensMan = (LensManager*)malloc(sizeof(LensManager));
    state->sourceMan = (SourceManager*)malloc(sizeof(SourceManager));
    InitLensManager(state->lensMan);
    InitSourceManager(state->sourceMan);

    for (int i=0; i<8; i++) {
        state->obsImg[i] = (SkyCoord){0.0f, 0.0f};
        state->obsImgIsVisible[i] = 0;
    }
    state->obsImgIsVisible[0] = 1; // By default, toggle the first image on

    state->images = (DynamicArrayImages*)malloc(sizeof(DynamicArrayImages));
    InitDAImages(state->images, 2);
    state->caustics = (DynamicArrayImages*)malloc(sizeof(DynamicArrayImages));
    InitDAImages(state->caustics, 1);
    state->criticalCurves = (DynamicArrayImages*)malloc(sizeof(DynamicArrayImages));
    InitDAImages(state->criticalCurves, 1);

    state->ctrSrcPlane = (PixelCoord){ctrSrcPlaneX, ctrSrcPlaneY};
    state->ctrImgPlane = (PixelCoord){ctrImgPlaneX, ctrImgPlaneY};

    state->grid_size = N_IGRIDP * N_IGRIDP;
    state->imgSGrid = (DynamicArraySkyCoord*)malloc(sizeof(DynamicArraySkyCoord));
    InitDynArraySkyCoord(state->imgSGrid, state->grid_size);
    state->srcSGrid = (DynamicArraySkyCoord*)malloc(sizeof(DynamicArraySkyCoord));
    InitDynArraySkyCoord(state->srcSGrid, state->grid_size);

    state->imgPlane = (ImagePlane*)malloc(sizeof(ImagePlane));
    InitImagePlane(state->imgPlane, img_fname, IMG_WIDTH_AS, IMG_HEIGHT_AS);
    state->gamma1 = 0;
    state->gamma2 = 0;

    TraceLog(LOG_DEBUG, "State initialized");
    return state;
}

void CleanStateGrids(State* state) {
    TraceLog(LOG_DEBUG, "Cleaning state grids with size %i", state->grid_size);
    CleanDynArraySkyCoord(state->imgSGrid, state->grid_size);
    CleanDynArraySkyCoord(state->srcSGrid, state->grid_size);
}

void FreeStateGrids(State* state) {
    CleanStateGrids(state);
    free(state->imgSGrid);
    free(state->srcSGrid);
}

void StateDestroy(State* state) {
    if (state != NULL) {
        FreeLensManager(state->lensMan);
        FreeSourceManager(state->sourceMan);

        FreeDAImages(state->images);
        free(state->images);
        FreeDAImages(state->caustics);
        free(state->caustics);
        FreeDAImages(state->criticalCurves);
        free(state->criticalCurves);

        //TODO: Free all the structures properly
        free(state->imgPlane);

        FreeStateGrids(state);
        free(state);
    }
}

void DisplayStateInfo(State* state) {
    // General info
    int zerx, zery, curx, cury;
    zerx = zeroExtraPanX + 15;
    zery = zeroExtraPanY + 15;
    curx = zerx; cury = zery;
    DrawText("Info:", curx, cury, 20, DARKMODE ? RAYWHITE : BLACKBOARD);
    // Lens info
    curx += 30; cury += 30;
    DrawText("Lens:", curx, cury, 20, DARKMODE ? RAYWHITE : BLACKBOARD);
    if (state->lensMan->count > 0) {
        Lens* lens = &state->lensMan->lenses[0];
        char lensposInfo[128];
        char lenspar1Info[128];
        char lenspar2Info[128];

        // TODO: Expand lens info based on model type
        curx += 10; cury += 25;
        snprintf(lensposInfo, sizeof(lensposInfo), "x = %.2f\", y = %.2f\"", lens->coord.RA, lens->coord.DEC);
        DrawText(lensposInfo, curx, cury, 16, DARKMODE ? RAYWHITE : BLACKBOARD);
        curx += 0; cury += 25;
        snprintf(lenspar1Info, sizeof(lenspar1Info), "b = %.2f\", qh = %.2f", lens->b, lens->qh);
        DrawText(lenspar1Info, curx, cury, 16, DARKMODE ? RAYWHITE : BLACKBOARD);
        curx += 0; cury += 25;
        snprintf(lenspar2Info, sizeof(lenspar2Info), "f = %.2f, PA = %.2f deg", lens->f, lens->PA * RAD2DEG);
        DrawText(lenspar2Info, curx, cury, 16, DARKMODE ? RAYWHITE : BLACKBOARD);
    } else {
        DrawText("No lens defined", zeroExtraPanX + 25, zeroExtraPanY + 55, 16, DARKMODE ? RAYWHITE : BLACKBOARD);
    }
    curx = zerx;;
    curx += 30; cury += 30;
    // Ext shear info
    DrawText("Ext Shear:", curx, cury, 20, DARKMODE ? RAYWHITE : BLACKBOARD);
    char extshearInfo[128];
    curx += 10; cury += 25;
    snprintf(extshearInfo, sizeof(extshearInfo), "g1 = %.2f, g2 = %.2f", state->gamma1, state->gamma2);
    DrawText(extshearInfo, curx, cury, 16, DARKMODE ? RAYWHITE : BLACKBOARD);

    curx = zerx; cury = zery;
    curx += 350; cury += 30;
    // Observed images info
    DrawText("Observed Images:", curx, cury, 20, DARKMODE ? RAYWHITE : BLACKBOARD);
    curx += 10;
    for (int i = 0; i < 8; i++) {
        if (state->obsImgIsVisible[i]) {
            char imgInfo[64];
            cury += 25;
            snprintf(imgInfo, sizeof(imgInfo), "Img %d: (%.2f\", %.2f\")", i + 1, state->obsImg[i].RA, state->obsImg[i].DEC);
            DrawText(imgInfo, curx, cury, 16, DARKMODE ? RAYWHITE : BLACKBOARD);
        }
    }
}

// **Einstein Radius**
//------------------------------------------------------------------------------------

double EinRad2PointMass(Lens* lens, AngCosmoDistances* cDists) {
    return 2 * (R_SCHWARTZ / MPC2METER) * lens->M * cDists->lens_dist_ratio / (ARCSEC2RADIAN * ARCSEC2RADIAN);
}

double EinRad2 (Lens* lens, Source* source) {
    // TODO: implement cosmological distances and remove prebaked DISTs
    // AngCosmoDistances* cDists = GetACosmoDistances(lens->zl, source->zs);

    AngCosmoDistances* cDists = (AngCosmoDistances*)malloc(sizeof(AngCosmoDistances));
    cDists->D_OL = DIST_OL;
    cDists->D_OS = DIST_OS;
    cDists->D_LS = DIST_LS;
    cDists->lens_dist_ratio = DIST_LS / (DIST_OL * DIST_OS);

    double result = 0.0;
    switch (lens->modelType) {
        case LTYPE_POINT:
            result = EinRad2PointMass(lens, cDists);
            break;

        case LTYPE_SIE:
            result = EinRad2PointMass(lens, cDists); //TODO: Implement SIE case
            break;

        case LTYPE_PERT:
            //TODO: Implement perturber lens case
            break;
        default:
            result = -1.;
            break;
    }
    free(cDists);
    return result;
}

// **Deflection Angles**
//------------------------------------------------------------------------------------

// **Point Lens**
//------------------------------------------------------------------------------------

// Function to calculate the deflection angle for a **point lens**
void CalcImgsPosPointLens(DynamicArrayImages* imagesArray, Lens* lens, Source* source) {
    double E2 = EinRad2(lens, source);
    double b1, b2, y, Q, xp1, xp2, xm1, xm2;
    double xSrcTmp, ySrcTmp;
    DynamicArraySkyCoord* imagesSCoord1 = malloc(sizeof(DynamicArraySkyCoord));
    DynamicArraySkyCoord* imagesSCoord2 = malloc(sizeof(DynamicArraySkyCoord));
    switch (source->modelType) {
        case STYPE_POINT:
            InitDynArraySkyCoord(imagesSCoord1, 1);
            InitDynArraySkyCoord(imagesSCoord2, 1);
            b1  = source->coord.RA  - lens->coord.RA;
            b2  = source->coord.DEC - lens->coord.DEC;
            y   = sqrt(b1*b1/E2+b2*b2/E2);
            Q   = sqrt(b1*b1/E2+b2*b2/E2 + 4) / y;
            xp1 = 0.5 * (1 + Q) * source->coord.RA;
            xp2 = 0.5 * (1 + Q) * source->coord.DEC;
            xm1 = 0.5 * (1 - Q) * source->coord.RA;
            xm2 = 0.5 * (1 - Q) * source->coord.DEC;
            AppendItemDynArraySkyCoord(imagesSCoord1, (SkyCoord){xp1, xp2});
            AppendItemDynArraySkyCoord(imagesSCoord2, (SkyCoord){xm1, xm2});
            AppendItemDAImages(imagesArray, *imagesSCoord1);
            AppendItemDAImages(imagesArray, *imagesSCoord2);
            break;
        case STYPE_EXT:
            InitDynArraySkyCoord(imagesSCoord1, N_EXTSRC);
            InitDynArraySkyCoord(imagesSCoord2, N_EXTSRC);
            for (int i = 0; i < N_EXTSRC; i++) {
                xSrcTmp = source->coord.RA  + source->radius * sin(angStep * i);
                ySrcTmp = source->coord.DEC + source->radius * cos(angStep * i);
                b1  = xSrcTmp - lens->coord.RA;
                b2  = ySrcTmp - lens->coord.DEC;
                y   = sqrt(b1*b1/E2+b2*b2/E2);
                Q   = sqrt(b1*b1/E2+b2*b2/E2 + 4) / y;
                xp1 = 0.5 * (1 + Q) * xSrcTmp;
                xp2 = 0.5 * (1 + Q) * ySrcTmp;
                xm1 = 0.5 * (1 - Q) * xSrcTmp;
                xm2 = 0.5 * (1 - Q) * ySrcTmp;
                AppendItemDynArraySkyCoord(imagesSCoord1, (SkyCoord){xp1, xp2});
                AppendItemDynArraySkyCoord(imagesSCoord2, (SkyCoord){xm1, xm2});
            }
            AppendItemDAImages(imagesArray, *imagesSCoord1);
            AppendItemDAImages(imagesArray, *imagesSCoord2);
            break;
        default:
            break;
    }
}

void tanCausticPoint(DynamicArrayImages* causticsArray, Lens* lens, Source* source){
    DynamicArraySkyCoord* tanSCoord = malloc(sizeof(DynamicArraySkyCoord));
    double phi;
    double Erd = sqrt(EinRad2(lens, source));
    InitDynArraySkyCoord(tanSCoord, N_EXTSRC);
    for (int i = 0; i < N_EXTSRC; i++) {
        phi = 2 * PI * ((double)i / N_EXTSRC);
        AppendItemDynArraySkyCoord(tanSCoord, (SkyCoord){lens->coord.RA  + Erd*cos(phi),
                                                         lens->coord.DEC + Erd*sin(phi)});
    }
    AppendItemDAImages(causticsArray, *tanSCoord);
}

void tanCriticalCurvePoint(DynamicArrayImages* criticalCurvesArray, Lens* lens, Source* source){
    DynamicArraySkyCoord* tanSCoord = malloc(sizeof(DynamicArraySkyCoord));
    double phi;
    double Erd = sqrt(EinRad2(lens, source));
    InitDynArraySkyCoord(tanSCoord, N_EXTSRC);
    for (int i = 0; i < N_EXTSRC; i++) {
        phi = 2 * PI * ((double)i / N_EXTSRC);
        AppendItemDynArraySkyCoord(tanSCoord, (SkyCoord){lens->coord.RA  + Erd*cos(phi),
                                                         lens->coord.DEC + Erd*sin(phi)});
    }
    AppendItemDAImages(criticalCurvesArray, *tanSCoord);
}

// **SIE Lens**
//------------------------------------------------------------------------------------


Vector2 Deflection_SIE(Lens *lens, double  xin, double  yin){
  double b = lens->b;
  double cospa = cos(lens->PA);
  double sinpa = sin(lens->PA);
  double dx = (xin-lens->coord.RA);
  double dy = (yin-lens->coord.DEC);

  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  dx*cospa + dy*sinpa;
  double y_t = -dx*sinpa + dy*cospa;

  if( fabs(x_t) < 0.0001 && fabs(y_t) < 0.0001 ){
    if( sign(x_t) ){
      x_t = -0.0001;
    } else {
      x_t =  0.0001;
    }
    if( sign(y_t) ){
      y_t = -0.0001;
    } else {
      y_t =  0.0001;
    }
  }

  double fac   = 1.0-lens->f*lens->f;
  double omega = lens->f*lens->f*x_t*x_t + y_t*y_t;
  // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  double fac2  = sqrt(fac/omega);

  double ax_t = (b/sqrt(fac))*atan (x_t*fac2);
  double ay_t = (b/sqrt(fac))*atanh(y_t*fac2);

  //rotate back according to position angle, no need to translate (this is equivalent to rotating by -pa using the same equations as above)
  Vector2 angle = {.x = ax_t*cospa - ay_t*sinpa, .y = ax_t*sinpa + ay_t*cospa};
  return angle;
}


Vector2 Deflection_Pert(Lens *lens, double  xin, double  yin){
    //TODO
    return (Vector2){0,0};
}


double PsiTilde(double phi, Lens* lens) {
    if (lens->f < 1.0) {
        double fp = sqrt(1 - lens->f*lens->f);
        double fr = sqrt(lens->f)/fp;
        double t1 = sin(phi - lens->PA) * asin( fp         * sin(phi - lens->PA));
        double t2 = cos(phi - lens->PA) * asinh(fp/lens->f * cos(phi - lens->PA));
        return fr * (t1 + t2);
    } else {
        return 1.0;
    }
}

double Psi(double x, double phi, Lens* lens) {
    return x * PsiTilde(phi, lens);
}

Vector2 AlphaSIE(double phi, Lens* lens) {
    if (lens->f < 1.0) {
        double fp = sqrt(1 - lens->f*lens->f);
        double fr = sqrt(lens->f)/fp;
        double a1 = fr * asinh(fp/lens->f * cos(phi));
        double a2 = fr * asin( fp         * sin(phi));
        return (Vector2){a1, a2};
    } else {
        return (Vector2){-cos(phi), -sin(phi)};
    }
}

double PhiFuncSIE(double phi, Vector2 yRot, Lens* lens){
    Vector2 a = AlphaSIE(phi, lens);
    return (a.x + yRot.x) * sin(phi) - (a.y + yRot.y) * cos(phi);
}

double radialDistImgSIE(Vector2 y, double phi, Lens* lens){
    return y.x*cos(phi) + y.y*sin(phi) + PsiTilde(phi+lens->PA, lens);;
}

// TODO: This should be a generic function that accepts func with more than one argument...
double RootFindBrentPhiFuncSIE(double x1, double x2, double tol, Vector2 yrot, Lens* lens) { //From Numerical Recipes - Sect 9.3
    const int ITMAX = 100;
    const double EPS = 1e-5;
    double a=x1,b=x2,c=x2,d,e,fa=PhiFuncSIE(a, yrot, lens),fb=PhiFuncSIE(b, yrot, lens),fc,p,q,r,s,tol1,xm;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        return -9999.;
    fc=fb;
    for (int iter=0;iter<ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) {
            return b;
        }
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=fabs(p);
            double min1=3.0*xm*q-fabs(tol1*q);
            double min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
            b += d;
        else {
            b += SIGN(tol1,xm);
            fb=PhiFuncSIE(b, yrot, lens);
        }
    }
    return -9999.;
}

void _CalcImgsPosSIELens(Vector2 y, Lens* lens, double Erd,
                         DynamicArraySkyCoord* imagesSCoord1,
                         DynamicArraySkyCoord* imagesSCoord2,
                         DynamicArraySkyCoord* imagesSCoord3,
                         DynamicArraySkyCoord* imagesSCoord4) {
    int N_PHI_ITER = 1000;
    double phi1, phi2, c1, c2, u, z, x;
    Vector2 yRot;
    SkyCoord s;
    // source position calculated from lens center
    yRot.x =  y.x - lens->coord.RA;
    yRot.y =  y.y - lens->coord.DEC;
    // source position in frame rotated to align with lens major axis
    yRot.x =  y.x * cos(lens->PA) + y.y * sin(lens->PA);
    yRot.y = -y.x * sin(lens->PA) + y.y * cos(lens->PA);

    // double y1_ =  y.x * cos(lens->PA) - y.y * sin(lens->PA);
    // double y2_ =  y.x * sin(lens->PA) + y.y * cos(lens->PA);
    int counter = 0;
    for (int i = 0; i < N_PHI_ITER - 1; i++) {
        phi1 = 2 * PI * ((double)i / N_PHI_ITER);
        phi2 = 2 * PI * ((double)(i+1) / N_PHI_ITER);
        c1 = PhiFuncSIE(phi1, yRot, lens);
        c2 = PhiFuncSIE(phi2, yRot, lens);
        if (sign(c1)+sign(c2)==0) {
            u = RootFindBrentPhiFuncSIE(phi1, phi2, 1e-5, yRot, lens);
            z = PhiFuncSIE(u, yRot, lens);
            if (z<0 && z>1e-5) continue;
            x = radialDistImgSIE(yRot, u, lens);
            if (x<=0) continue;
            s = MapSkyPolarToSkyCartesian(x * Erd, u + lens->PA);
#if 0
//External shear
            Vector2 tmpalpha, tmp_x;
            tmpalpha.x = s.RA - y.x;
            tmpalpha.y = s.DEC - y.y;
            tmp_x.x = s.RA;
            tmp_x.y = s.DEC;
            s.RA  = (1.0-state->gamma1)*y.x - state->gamma2*y.y - tmpalpha.x;
            s.DEC = (1.0+state->gamma1)*y.y - state->gamma2*y.x - tmpalpha.y;
#endif
            switch (counter) {
                case 0:
                    AppendItemDynArraySkyCoord(imagesSCoord1, s);
                    break;
                case 1:
                    AppendItemDynArraySkyCoord(imagesSCoord2, s);
                    break;
                case 2:
                    AppendItemDynArraySkyCoord(imagesSCoord3, s);
                    break;
                case 3:
                    AppendItemDynArraySkyCoord(imagesSCoord4, s);
                    break;
                default:
                    break;
                }
                counter++;
        }
    }
}

void CalcImgsPosSIELens(DynamicArrayImages* imagesArray, Lens* lens, Source* source) {
    double Erd = sqrt(EinRad2(lens, source));
    Vector2 y;
    DynamicArraySkyCoord* imagesSCoord1 = malloc(sizeof(DynamicArraySkyCoord));
    DynamicArraySkyCoord* imagesSCoord2 = malloc(sizeof(DynamicArraySkyCoord));
    DynamicArraySkyCoord* imagesSCoord3 = malloc(sizeof(DynamicArraySkyCoord));
    DynamicArraySkyCoord* imagesSCoord4 = malloc(sizeof(DynamicArraySkyCoord));
    TraceLog(LOG_DEBUG, "Calculating images for SIE lens\n");
    switch (source->modelType) {
        case STYPE_POINT:
            InitDynArraySkyCoord(imagesSCoord1, 1);
            InitDynArraySkyCoord(imagesSCoord2, 1);
            InitDynArraySkyCoord(imagesSCoord3, 1);
            InitDynArraySkyCoord(imagesSCoord4, 1);
            TraceLog(LOG_DEBUG, "   Init DAimages for SIE lens\n");
            y = (Vector2){source->coord.RA/Erd, source->coord.DEC/Erd};
            // source position calculated from lens center
            y.x =  y.x - lens->coord.RA;
            y.y =  y.y - lens->coord.DEC;
            _CalcImgsPosSIELens(y, lens, Erd, imagesSCoord1, imagesSCoord2, imagesSCoord3, imagesSCoord4);
            TraceLog(LOG_DEBUG, "   Appending images to array\n");
            AppendItemDAImages(imagesArray, *imagesSCoord1);
            AppendItemDAImages(imagesArray, *imagesSCoord2);
            AppendItemDAImages(imagesArray, *imagesSCoord3);
            AppendItemDAImages(imagesArray, *imagesSCoord4);
            break;
        case STYPE_EXT:
            InitDynArraySkyCoord(imagesSCoord1, N_EXTSRC);
            InitDynArraySkyCoord(imagesSCoord2, N_EXTSRC);
            InitDynArraySkyCoord(imagesSCoord3, N_EXTSRC);
            InitDynArraySkyCoord(imagesSCoord4, N_EXTSRC);
            TraceLog(LOG_DEBUG, "   Init DAimages for SIE lens\n");
            Vector2 angle; // Debugging variable
            SkyCoord s; // Debugging variable
            for (int j = 0; j < N_EXTSRC; j++) {
                // source position calculated from lens center
                y = (Vector2){source->coord.RA/Erd  + source->radius * sin(angStep * j),
                source->coord.DEC/Erd + source->radius * cos(angStep * j)};
#if 0
                TraceLog(LOG_DEBUG, "   Calculating deflection\n");
                angle = defl_SIE(Erd, lens, y.x, y.y);
                TraceLog(LOG_DEBUG, "   Calculated deflection\n");
                s.RA = y.x - angle.x;
                s.DEC = y.y - angle.y;
                s.RA  = s.RA  - lens->coord.RA;
                s.DEC = s.DEC - lens->coord.DEC;
                AppendItemDynArraySkyCoord(imagesSCoord1, s);
                AppendItemDynArraySkyCoord(imagesSCoord2, s);
                AppendItemDynArraySkyCoord(imagesSCoord3, s);
                AppendItemDynArraySkyCoord(imagesSCoord4, s);
            }
            TraceLog(LOG_DEBUG, "   Appending images to array\n");
            TraceLog(LOG_DEBUG, "last angle (%d, %d)", angle.x, angle.y);
            TraceLog(LOG_DEBUG, "last SkyCoord (%d, %d)", s.RA, s.DEC);
            TraceLog(LOG_ERROR, "---->HERE:  Length of imagesSCoord1 is %d\n", imagesSCoord1->size);
#else
                y.x =  y.x - lens->coord.RA;
                y.y =  y.y - lens->coord.DEC;
                _CalcImgsPosSIELens(y, lens, Erd, imagesSCoord1, imagesSCoord2, imagesSCoord3, imagesSCoord4);
            }
#endif
            AppendItemDAImages(imagesArray, *imagesSCoord1);
            AppendItemDAImages(imagesArray, *imagesSCoord2);
            AppendItemDAImages(imagesArray, *imagesSCoord3);
            AppendItemDAImages(imagesArray, *imagesSCoord4);
            break;
        default:
            break;
    }
}

void cutSIE(DynamicArrayImages* causticsArray, Lens* lens, Source* source){
    DynamicArraySkyCoord* cutSCoord = malloc(sizeof(DynamicArraySkyCoord));
    Vector2 a;
    SkyCoord y;
    double phi;
    double Erd = sqrt(EinRad2(lens, source));
    InitDynArraySkyCoord(cutSCoord, N_EXTSRC);
    for (int i = 0; i < N_EXTSRC; i++) {
        phi = 2 * PI * ((double)i / N_EXTSRC);
        a = AlphaSIE(phi, lens);
        y.RA  = a.x * cos(lens->PA) - a.y * sin(lens->PA);
        y.DEC = a.x * sin(lens->PA) + a.y * cos(lens->PA);
        // centred on the lens
        y.RA  =  -(y.RA*Erd  + lens->coord.RA);
        y.DEC =  -(y.DEC*Erd + lens->coord.DEC);
        AppendItemDynArraySkyCoord(cutSCoord, y);
    }
    AppendItemDAImages(causticsArray, *cutSCoord);
}

void tanCausticSIE(DynamicArrayImages* causticsArray, Lens* lens, Source* source){
    DynamicArraySkyCoord* tanSCoord = malloc(sizeof(DynamicArraySkyCoord));
    SkyCoord y;
    Vector2 a;
    double phi, delta, y1_, y2_;
    double Erd = sqrt(EinRad2(lens, source));
    double sqrtf = sqrt(lens->f);
    InitDynArraySkyCoord(tanSCoord, N_EXTSRC);
    for (int i = 0; i < N_EXTSRC; i++) {
        phi = 2 * PI * ((double)i / N_EXTSRC);
        delta = sqrt(cos(phi)*cos(phi)+lens->f*lens->f*sin(phi)*sin(phi));
        a = AlphaSIE(phi, lens);
        y1_ = sqrtf/delta*cos(phi)-a.x;
        y2_ = sqrtf/delta*sin(phi)-a.y;
        y.RA   = y1_ * cos(lens->PA) - y2_ * sin(lens->PA);
        y.DEC  = y1_ * sin(lens->PA) + y2_ * cos(lens->PA);
        // centred on the lens
        y.RA  =  y.RA *Erd - lens->coord.RA;
        y.DEC =  y.DEC*Erd - lens->coord.DEC;
        AppendItemDynArraySkyCoord(tanSCoord, y);
    }
    AppendItemDAImages(causticsArray, *tanSCoord);
}

void tanCriticalCurveSIE(DynamicArrayImages* criticalCurveArray, Lens* lens, Source* source){
    DynamicArraySkyCoord* tanSCoord = malloc(sizeof(DynamicArraySkyCoord));
    SkyCoord x;
    Vector2 a;
    double phi, delta, r, x1, x2;
    double Erd = sqrt(EinRad2(lens, source));
    InitDynArraySkyCoord(tanSCoord, N_EXTSRC);
    for (int i = 0; i < N_EXTSRC; i++) {
        phi = 2 * PI * ((double)i / N_EXTSRC);
        delta = sqrt(cos(phi)*cos(phi)+lens->f*lens->f*sin(phi)*sin(phi));
        r  = sqrt(lens->f)/delta;
        x.RA  = r * cos(phi + lens->PA);
        x.DEC = r * sin(phi + lens->PA);
        // centred on the lens
        x.RA  =  x.RA *Erd - lens->coord.RA;
        x.DEC =  x.DEC*Erd - lens->coord.DEC;
        AppendItemDynArraySkyCoord(tanSCoord, x);
    }
    AppendItemDAImages(criticalCurveArray, *tanSCoord);
}

// **SPEMD Lens**
//------------------------------------------------------------------------------------

Vector2 Deflection_SPEMD(Lens *lens, double x_in, double y_in) {
    double b, rc2;
	double xtr, ytr;
	double dfl[2];
    double w, w2, w2p1q, q, gamma, pgam, A, B, Re;
    double alpha_x = 0.0;
	double alpha_y = 0.0;

    // translate and rotate into the frame of the profile
	physics_util_translate_rotate_fwd(x_in, y_in, (double)lens->coord.RA, (double)lens->coord.DEC,
			lens->sin_PA, lens->cos_PA, &xtr, &ytr);

    // b defined as theta_E, see Enzi et al 2020
    Re = lens->b;
    q = lens->f;
    w = lens->rc / Re;
    gamma = 2*lens->qh+1;
    w2 = w*w;
    w2p1q = w*w + 1/q;
    pgam = (3-gamma)/2 ;
    A = (1.0/2 * (4 - gamma)/(3 - gamma) * pow(q,1.0/2));
    B = (pow(w2p1q, pgam) - pow(w2, pgam));
    b = pow(Re, gamma - 1 ) / A / B;
    b =  ( 1.5 - lens->qh) * (b/(2.0*sqrt(lens->f)));

    // TODO: need to re-scale based on the critical density,
	// since b is the Einstein radius, which has Ds and Dls baked in

    // call fastell with the correct units
	rc2 = lens->rc * lens->rc;
	fastelldefl_(&xtr, &ytr, &b, &lens->qh, &lens->f, &rc2, dfl);

	// rotate the deflection angles back
	physics_util_rotate_bwd(dfl[0], dfl[1], lens->sin_PA, lens->cos_PA, &alpha_x, &alpha_y);
    return (Vector2){alpha_x, alpha_y};
}


// **Generic Deflection**
//------------------------------------------------------------------------------------


Vector2 Deflection(Lens *lens, double  xin, double  yin){
    switch(lens->modelType){
      case LTYPE_POINT:
        return (Vector2){0,0};
        break;
      case LTYPE_SIE:
        return Deflection_SIE(lens, xin, yin);
        break;
      case LTYPE_SPEMD:
        return Deflection_SPEMD(lens, xin, yin);
        break;
      case LTYPE_PERT:
        return Deflection_Pert(lens, xin, yin);
        break;
      default:
        return (Vector2){0,0};
        break;
    }
}

SkyCoord SingleSourcePosition(State *state, Lens *lens, SkyCoord imgPos){
    double src_x = 0.0;
	double src_y = 0.0;
    Vector2 angle = Deflection(lens, imgPos.RA, imgPos.DEC);
    src_x = (1.0-state->gamma1)*imgPos.RA - state->gamma2*imgPos.DEC - (double)angle.x;
    src_y = (1.0+state->gamma1)*imgPos.DEC - state->gamma2*imgPos.RA - (double)angle.y;
    return (SkyCoord){(float)src_x, (float)src_y};
}

void Deflection_Total_OnAllPixels(State *state){
    double dumx = 0.0;
    double dumy = 0.0;
    double xin,yin;

    Vector2 angle;
    double ax, ay;
    for(int j=0;j<state->imgPlane->Nm;j++){
        xin = state->imgPlane->x[j];
        yin = state->imgPlane->y[j];
        ax   = 0.0;
        ay   = 0.0;
        for(int i=0;i<state->lensMan->count;i++){
            angle = Deflection(&state->lensMan->lenses[i], xin, yin);
            ax += angle.x;
            ay += angle.y;
        }
        state->imgPlane->defl_x[j] = (1.0-state->gamma1)*xin - state->gamma2*yin - ax;
        state->imgPlane->defl_y[j] = (1.0+state->gamma1)*yin - state->gamma2*xin - ay;
    }
}


// **Image Positions**
//------------------------------------------------------------------------------------

void CalcImgsPos (DynamicArrayImages* imagesArray, Lens* lens, Source* source) {
    // TraceLog(LOG_DEBUG, "Calculating images positions");
    switch (lens->modelType) {
        case LTYPE_POINT:

            CalcImgsPosPointLens(imagesArray, lens, source);
            break;

        case LTYPE_SIE:
            TraceLog(LOG_DEBUG, "------- SIE img pos --------");
            CalcImgsPosSIELens(imagesArray, lens, source);
            break;

        case LTYPE_SPEMD:
            TraceLog(LOG_DEBUG, "------- SPEMD img pos --------");
            // CalcImgsPosSPEMDLens(imagesArray, lens, source);
            break;

        case LTYPE_PERT:
            //TODO: Implement perturber lens case
            break;

        default:
            break;
    }
}

// **Critical Lines and Caustics**
//------------------------------------------------------------------------------------

void CalcCaustics (DynamicArrayImages* causticsArray, Lens* lens, Source* source) {
    TraceLog(LOG_DEBUG, "Calculating caustics");
    switch (lens->modelType) {
        case LTYPE_POINT:
            tanCausticPoint(causticsArray, lens, source);
            break;

        case LTYPE_SIE:
            tanCausticSIE(causticsArray, lens, source);
            cutSIE(causticsArray, lens, source);
            break;

        case LTYPE_SPEMD:
            //TODO:
            // tanCausticSPEMD(causticsArray, lens, source);
            // cutSPEMD(causticsArray, lens, source);
            break;

        case LTYPE_PERT:
            //TODO: Implement perturber lens case
            break;

        default:
            break;
    }
}

void CalcCriticalCurves (DynamicArrayImages* CriticalCurvesArray, Lens* lens, Source* source) {
    // TraceLog(LOG_DEBUG, "Calculating critical curves");
    switch (lens->modelType) {
        case LTYPE_POINT:
            tanCriticalCurvePoint(CriticalCurvesArray, lens, source);
            // TraceLog(LOG_DEBUG, "Calculated critical curves for point lens");
            break;

        case LTYPE_SIE:
            tanCriticalCurveSIE(CriticalCurvesArray, lens, source);
            // TraceLog(LOG_DEBUG, "Calculated critical curves for SIE lens");
            break;

        case LTYPE_SPEMD:
            // TODO:
            // tanCriticalCurveSPEMD(CriticalCurvesArray, lens, source);
            // TraceLog(LOG_DEBUG, "Calculated critical curves for SPEMD lens");
            break;

        case LTYPE_PERT:
            //TODO: Implement perturber lens case
            break;

        default:
            break;
    }
}


// **Draw functions**
//------------------------------------------------------------------------------------

void DrawSource(Source* source) {
    PixelCoord p1, p2;
    double xtmp, ytmp;
    switch (source->modelType) {
        case STYPE_POINT:
            p1 = MapSkyToPixel(source->coord, screenWidth, screenHeight, 1);
            if (CheckPixelInTargetWindow(p1, 1)){
                DrawCircle(p1.x, p1.y, 2, ORANGE);
            }
            break;
        case STYPE_EXT:
            for (int i = 0; i < N_EXTSRC - 1; i++) {
                xtmp = source->coord.RA  + source->radius * sin(angStep * i);
                ytmp = source->coord.DEC + source->radius * cos(angStep * i);
                p1 = MapSkyToPixel((SkyCoord){xtmp, ytmp}, screenWidth, screenHeight, 1);
                xtmp = source->coord.RA  + source->radius * sin(angStep * (i + 1));
                ytmp = source->coord.DEC + source->radius * cos(angStep * (i + 1));
                p2 = MapSkyToPixel((SkyCoord){xtmp, ytmp}, screenWidth, screenHeight, 1);
                if (CheckPixelInTargetWindow(p1, 1) && CheckPixelInTargetWindow(p2, 1)){
                    if (DRAW_EXT_AS_POINTS) {
                        DrawCircle(p1.x, p1.y, 1.2, ORANGE);
                        DrawCircle(p2.x, p2.y, 1.2, ORANGE);
                    } else {
                        DrawLineEx((Vector2){p1.x, p1.y}, (Vector2){p2.x, p2.y}, 1.2, ORANGE);
                    }
                }
            }
            xtmp = source->coord.RA  + source->radius * sin(angStep * N_EXTSRC);
            ytmp = source->coord.DEC + source->radius * cos(angStep * N_EXTSRC);
            p1 = MapSkyToPixel((SkyCoord){xtmp, ytmp}, screenWidth, screenHeight, 1);
            xtmp = source->coord.RA; // + source->radius * 0;
            ytmp = source->coord.DEC + source->radius;
            p2 = MapSkyToPixel((SkyCoord){xtmp, ytmp}, screenWidth, screenHeight, 1);
            if (CheckPixelInTargetWindow(p1, 1) && CheckPixelInTargetWindow(p2, 1)){
                if (DRAW_EXT_AS_POINTS) {
                    DrawCircle(p1.x, p1.y, 1.2, ORANGE);
                    DrawCircle(p2.x, p2.y, 1.2, ORANGE);
                } else {
                    DrawLineEx((Vector2){p1.x, p1.y}, (Vector2){p2.x, p2.y}, 1.2, ORANGE);
                }
            }
            break;
        default:
            break;
    }
}

void DrawAllSources(State* state){
    for (int i = 0; i < state->sourceMan->count; i++) {
        DrawSource(&state->sourceMan->sources[i]);
    }
}

void DrawPointImages(SkyCoord SCoord, int colorIndex, float sizePoint, bool IsSourcePlane){
    PixelCoord p1;
    p1 = MapSkyToPixel(SCoord, screenWidth, screenHeight, IsSourcePlane);
    if (CheckPixelInTargetWindow(p1, IsSourcePlane)){
        DrawCircle(p1.x, p1.y, sizePoint, colors[colorIndex % numColors]);
    }
}

void DrawImages(DynamicArraySkyCoord imgsSCoordArray, int colorIndex){
    PixelCoord p1, p2;
    float sizePoint = 2.0;
    if (!(imgsSCoordArray.size<=0)) {
        if (imgsSCoordArray.size==1) {
            p1 = MapSkyToPixel(imgsSCoordArray.coords[0], screenWidth, screenHeight, 0);
            if (CheckPixelInTargetWindow(p1, 0)){
                DrawCircle(p1.x, p1.y, sizePoint, colors[colorIndex % numColors]);
            }
        } else {
            for (int i = 0; i < imgsSCoordArray.size - 1; i++) {
                p1 = MapSkyToPixel(imgsSCoordArray.coords[i], screenWidth, screenHeight, 0);
                p2 = MapSkyToPixel(imgsSCoordArray.coords[i + 1], screenWidth, screenHeight, 0);
                if (CheckPixelInTargetWindow(p1, 0) && CheckPixelInTargetWindow(p2, 0)){
                    if (DRAW_EXT_AS_POINTS) {
                        DrawCircle(p1.x, p1.y, sizePoint*0.6, ORANGE);
                        DrawCircle(p2.x, p2.y, sizePoint*0.6, ORANGE);
                    } else {
                        DrawLineEx((Vector2){p1.x, p1.y}, (Vector2){p2.x, p2.y}, 1.2, colors[colorIndex % numColors]);
                    }
                }
            }
            p1 = MapSkyToPixel(imgsSCoordArray.coords[imgsSCoordArray.size - 1], screenWidth, screenHeight, 0);
            p2 = MapSkyToPixel(imgsSCoordArray.coords[0], screenWidth, screenHeight, 0);
            if (CheckPixelInTargetWindow(p1, 0) && CheckPixelInTargetWindow(p2, 0)){
                if (DRAW_EXT_AS_POINTS) {
                    DrawCircle(p1.x, p1.y, sizePoint*0.6, ORANGE);
                    DrawCircle(p2.x, p2.y, sizePoint*0.6, ORANGE);
                } else {
                    DrawLineEx((Vector2){p1.x, p1.y}, (Vector2){p2.x, p2.y}, sizePoint*0.6, colors[colorIndex % numColors]);
                }
            }
        }
    }
}

void DrawAllImages(State* state){
    if (state->images->size > 0){
        for (int i = 0; i < state->images->size; i++) {
            DrawImages(state->images->imgSCoordDArray[i], i);
        }
    }
}

void DrawAllObsImages(State* state){
    Vector2 angle;
    SkyCoord s;
    for (int i = 0; i < 8; i++) {
        if(state->obsImgIsVisible[i]) {
            DrawPointImages(state->obsImg[i], i, 5.0, false);
            s = SingleSourcePosition(state, &state->lensMan->lenses[0], state->obsImg[i]);
            DrawPointImages(s, i, 5.0, true);
            // TraceLog(LOG_DEBUG, "SingleSourcePosition: imgPos (%.5f, %.5f) -> srcPos (%.5f, %.5f)", state->obsImg[i].RA, state->obsImg[i].DEC, s.RA, s.DEC);
        }
    }
}

void DrawCenterPaLens(State* state){
    SkyCoord scenter;
    SkyCoord s1, s2, s3, s4;

    PixelCoord pcenter;
    PixelCoord p1, p2, p3, p4;

    double crosshairSize = 0.15; // arcsec
    double f; // axis ratio
    Vector2 dir;

    for (int i = 0; i < state->lensMan->count; i++) {
        scenter = state->lensMan->lenses[i].coord;
        f = state->lensMan->lenses[i].f;

        pcenter = MapSkyToPixel(scenter, screenWidth, screenHeight, 0);
        dir.x = state->lensMan->lenses[i].cos_PA;
        dir.y = state->lensMan->lenses[i].sin_PA;

        s1 = (SkyCoord){scenter.RA - crosshairSize * dir.y, scenter.DEC + crosshairSize * dir.x};
        s2 = (SkyCoord){scenter.RA + crosshairSize * dir.y, scenter.DEC - crosshairSize * dir.x};
        s3 = (SkyCoord){scenter.RA - f * crosshairSize * dir.x, scenter.DEC - f * crosshairSize * dir.y};
        s4 = (SkyCoord){scenter.RA + f * crosshairSize * dir.x, scenter.DEC + f * crosshairSize * dir.y};
        p1 = MapSkyToPixel(s1, screenWidth, screenHeight, 0);
        p2 = MapSkyToPixel(s2, screenWidth, screenHeight, 0);
        p3 = MapSkyToPixel(s3, screenWidth, screenHeight, 0);
        p4 = MapSkyToPixel(s4, screenWidth, screenHeight, 0);
        DrawLineEx((Vector2){p1.x, p1.y}, (Vector2){p2.x, p2.y}, 1.2, DARKMODE ? WHITE : BLACK);
        DrawLineEx((Vector2){p3.x, p3.y}, (Vector2){p4.x, p4.y}, 1.2, DARKMODE ? WHITE : BLACK);
    }
}

void DrawCaustics(DynamicArraySkyCoord causticSCoordArray, int colorIndex) {
    PixelCoord p1, p2;
    if (causticSCoordArray.size > 1) {
        for (int i = 0; i < causticSCoordArray.size - 1; i++) {
            p1 = MapSkyToPixel(causticSCoordArray.coords[i], screenWidth, screenHeight, 1);
            p2 = MapSkyToPixel(causticSCoordArray.coords[i + 1], screenWidth, screenHeight, 1);
            if (CheckPixelInTargetWindow(p1, 1) && CheckPixelInTargetWindow(p2, 1)){
                DrawLineEx((Vector2){p1.x, p1.y}, (Vector2){p2.x, p2.y}, 1.2, colors[colorIndex % numColors]);
            }
        }
        p1 = MapSkyToPixel(causticSCoordArray.coords[causticSCoordArray.size - 1], screenWidth, screenHeight, 1);
        p2 = MapSkyToPixel(causticSCoordArray.coords[0], screenWidth, screenHeight, 1);
        if (CheckPixelInTargetWindow(p1, 1) && CheckPixelInTargetWindow(p2, 1)){
            DrawLineEx((Vector2){p1.x, p1.y}, (Vector2){p2.x, p2.y}, 1.2, colors[colorIndex % numColors]);
        }
    }
}

void DrawCriticalCurves(DynamicArraySkyCoord criticalCurvesSCoordArray, int colorIndex) {
    PixelCoord p1, p2;
    if (criticalCurvesSCoordArray.size==1) {
        p1 = MapSkyToPixel(criticalCurvesSCoordArray.coords[0], screenWidth, screenHeight, 0);
        if (CheckPixelInTargetWindow(p1, 0)){
            DrawCircle(p1.x, p1.y, 2, colors[colorIndex % numColors]);
        }
    } else {
        for (int i = 0; i < criticalCurvesSCoordArray.size - 1; i++) {
            p1 = MapSkyToPixel(criticalCurvesSCoordArray.coords[i], screenWidth, screenHeight, 0);
            p2 = MapSkyToPixel(criticalCurvesSCoordArray.coords[i + 1], screenWidth, screenHeight, 0);
            if (CheckPixelInTargetWindow(p1, 0) && CheckPixelInTargetWindow(p2, 0)){
                DrawLineEx((Vector2){p1.x, p1.y}, (Vector2){p2.x, p2.y}, 1.2, LIME);
            }
        }
        p1 = MapSkyToPixel(criticalCurvesSCoordArray.coords[criticalCurvesSCoordArray.size - 1], screenWidth, screenHeight, 0);
        p2 = MapSkyToPixel(criticalCurvesSCoordArray.coords[0], screenWidth, screenHeight, 0);
        if (CheckPixelInTargetWindow(p1, 0) && CheckPixelInTargetWindow(p2, 0)){
            DrawLineEx((Vector2){p1.x, p1.y}, (Vector2){p2.x, p2.y}, 1.2, LIME);
        }
    }
}

void DrawAllCausticAndCriticalLine(State* state){
    for (int i = 0; i < state->caustics->size; i++)       DrawCaustics(state->caustics->imgSCoordDArray[i], i);
    for (int i = 0; i < state->criticalCurves->size; i++) DrawCriticalCurves(state->criticalCurves->imgSCoordDArray[i], i);
}

// **Grids on Source and Image Planes**
//------------------------------------------------------------------------------------

void CalcImgPlaneGrid(DynamicArraySkyCoord* imgSGrid, Lens* lens, Source* source) {
    TraceLog(LOG_DEBUG, "Setting img plane grid of size %i", imgSGrid->size);
    for (int i = 0; i < N_IGRIDP; i++) {
        for (int j = 0; j < N_IGRIDP; j++) {
            double RAi  = (double)((double)(i-N_IGRIDP/2)/N_IGRIDP) * 10;
            double DECi = (double)((double)(j-N_IGRIDP/2)/N_IGRIDP) * 10;
            AppendItemDynArraySkyCoord(imgSGrid, (SkyCoord){RAi, DECi});
        }
    }
}

void CalcSrcPlaneGrid(DynamicArraySkyCoord* srcSGrid, State* state, Lens* lens) {
    TraceLog(LOG_DEBUG, "Setting src plane grid of size %i", srcSGrid->size);
    SkyCoord y, yRot;
    Vector2 a;
    double x1, x2, x1Rot, x2Rot, eRad2, eRad;
    eRad2 = EinRad2(lens, &state->sourceMan->sources[0]);
    eRad  = sqrt(eRad2);
    for (int i = 0; i < state->grid_size; i++) {
        x1 = state->imgSGrid->coords[i].RA;
        x2 = state->imgSGrid->coords[i].DEC;
        x1Rot = x1 * cos(-lens->PA) - x2 * sin(-lens->PA); //TODO: Use a 2d rotation function (also in case LTYPE_SIE)
        x2Rot = x1 * sin(-lens->PA) + x2 * cos(-lens->PA);
        switch(lens->modelType) {
            case LTYPE_POINT: //TODO: This is wrong (?)
                y.RA  = x1 - eRad2 / x1;
                y.DEC = x2 - eRad2 / x2;
                break;
            case LTYPE_SIE:
                a = AlphaSIE(atan2(x2Rot, x1Rot), lens);
#if 1
                //TODO: Check order
                yRot.RA  = (1.0-state->gamma1)*x1Rot - state->gamma2*x2Rot - a.x * eRad;
                yRot.DEC = (1.0+state->gamma1)*x2Rot - state->gamma2*x1Rot - a.y * eRad;
                y.RA  = yRot.RA * cos(lens->PA) - yRot.DEC * sin(lens->PA);
                y.DEC = yRot.RA * sin(lens->PA) + yRot.DEC * cos(lens->PA);
#else
                yRot.RA  = x1Rot - a.x * eRad;
                yRot.DEC = x2Rot - a.y * eRad;
                y.RA  = yRot.RA * cos(lens->PA) - yRot.DEC * sin(lens->PA);
                y.DEC = yRot.RA * sin(lens->PA) + yRot.DEC * cos(lens->PA);
#endif
                break;
            case LTYPE_PERT:
                //TODO: Implement perturber lens case
                break;
            default:
                break;
        }
        AppendItemDynArraySkyCoord(srcSGrid, y);
    }
}

void DrawDebugGridColoredPointAtPixelCoord(PixelCoord pixel, DynamicArraySkyCoord* skyGrid, int i, double gridLimDistance) {
    float distance = sqrt(skyGrid->coords[i].RA * skyGrid->coords[i].RA + skyGrid->coords[i].DEC * skyGrid->coords[i].DEC);
    float angle = atan2(skyGrid->coords[i].DEC, skyGrid->coords[i].RA);
    float maxDistance = 10.0; // Assuming the maximum distance is 10 arcsecs
    float hue = fmod((distance / maxDistance) * 360.0 * 3, 360.0); // Map distance to hue
    // hue = hue * cos(angle) * sin(angle); // Map distance and angle to hue
    Color pointColor = ColorFromHSV(hue, 1.0, 1.0); // Full saturation and value for rainbow colors
    if (distance < gridLimDistance) DrawCircle(pixel.x, pixel.y, 1.2, pointColor);
}

void DrawSkyGrid(DynamicArraySkyCoord* sGrid, bool IsSrcPlane, double gridLimDistance) {
    if (sGrid == NULL || sGrid->coords == NULL) return;
    for (int i = 0; i < sGrid->size; i++) {
        PixelCoord _imgPGrid = MapSkyToPixel(sGrid->coords[i], screenWidth, screenHeight, IsSrcPlane);
        if (CheckPixelInTargetWindow(_imgPGrid, IsSrcPlane)){
            DrawDebugGridColoredPointAtPixelCoord(_imgPGrid, sGrid, i, gridLimDistance);
        }
    }
}
void DrawSkyGridSrc(DynamicArraySkyCoord* sGrid, DynamicArraySkyCoord* iGrid, double gridLimDistance) {
    if (sGrid == NULL || sGrid->coords == NULL) return;
    for (int i = 0; i < iGrid->size; i++) {
        PixelCoord _srcPGrid = MapSkyToPixel(sGrid->coords[i], screenWidth, screenHeight, 1);
        if (CheckPixelInTargetWindow(_srcPGrid, 1)){
            DrawDebugGridColoredPointAtPixelCoord(_srcPGrid, iGrid, i, gridLimDistance);
        }
    }
}

// **Updating Models**
void UpdateLensAndSourceModels(State* state){
    // Mode 1 : produce images from source plane positions
    CleanDAImages(state->images, 2);
    TraceLog(LOG_DEBUG, "------- calc img pos --------");
    CalcImgsPos(state->images, &state->lensMan->lenses[0], &state->sourceMan->sources[0]);
    TraceLog(LOG_DEBUG, "---- [done] calc img pos ----");

    CleanDAImages(state->caustics, 1);
    CalcCaustics(state->caustics, &state->lensMan->lenses[0], &state->sourceMan->sources[0]);
    CleanDAImages(state->criticalCurves, 1);
    CalcCriticalCurves(state->criticalCurves, &state->lensMan->lenses[0], &state->sourceMan->sources[0]);

    // Mode 2 : make SIE lens with random mass and f
    CleanStateGrids(state);
    CalcImgPlaneGrid(state->imgSGrid, &state->lensMan->lenses[0], &state->sourceMan->sources[0]);
    CalcSrcPlaneGrid(state->srcSGrid, state, &state->lensMan->lenses[0]);
}

// **GUI**
//------------------------------------------------------------------------------------
typedef struct {
    bool formActive;

    LensModelType lensModelType;
    float lensMassInput;
    float lensbInput;
    float lensqhInput;
    float lensrcInput;
    float lensFInput;
    float lensPAInput;
    float gamma_1;
    float gamma_2;


    SourceModelType sourceModelType;
    float sourceRadius;

    SkyCoord lensPositionSlider;
    SkyCoord sourcePositionSlider;

    SkyCoord obsImg1PositionSlider;
    SkyCoord obsImg2PositionSlider;
    SkyCoord obsImg3PositionSlider;
    SkyCoord obsImg4PositionSlider;
    SkyCoord obsImg5PositionSlider;
    SkyCoord obsImg6PositionSlider;
    SkyCoord obsImg7PositionSlider;
    SkyCoord obsImg8PositionSlider;

    bool obsImg1IsVisible;
    bool obsImg2IsVisible;
    bool obsImg3IsVisible;
    bool obsImg4IsVisible;
    bool obsImg5IsVisible;
    bool obsImg6IsVisible;
    bool obsImg7IsVisible;
    bool obsImg8IsVisible;

    float tmpSrcScale, tmpImgScale;
    bool  viewVorPoints, viewFITS, viewScaleBar, viewParity, viewExtSrcPoints;
    float gridLimDistance;
} FormGUI; //TODO: I should probably put this inside State

void UpdateLensAndSourceModelsFromGUI(State* state, FormGUI* formGUI) {
    if (state->lensMan->count > 0) {
        Lens* lens = &state->lensMan->lenses[0];
        lens->modelType = formGUI->lensModelType;
        SELECTED_LENS_MODEL = dropdownLensModel;
        lens->M  = formGUI->lensMassInput * 1e12;
        lens->b  = (double)formGUI->lensbInput;
        lens->qh = (double)formGUI->lensqhInput;
        lens->rc = (double)formGUI->lensrcInput;
        lens->f  = formGUI->lensFInput;
        lens->PA = formGUI->lensPAInput;
        lens->cos_PA = cos(lens->PA);
        lens->sin_PA = sin(lens->PA);
        lens->coord.RA = formGUI->lensPositionSlider.RA;
        lens->coord.DEC = formGUI->lensPositionSlider.DEC;
    }
    state->gamma1 = formGUI->gamma_1;
    state->gamma2 = formGUI->gamma_2;

    if (state->sourceMan->count > 0) {
        Source* source = &state->sourceMan->sources[0];
        source->modelType = formGUI->sourceModelType;
        source->coord.RA = formGUI->sourcePositionSlider.RA;
        source->coord.DEC = formGUI->sourcePositionSlider.DEC;
        if (formGUI->sourceModelType == STYPE_EXT) {
            source->radius = formGUI->sourceRadius;
        }
    }

    state->obsImg[1-1].RA  = formGUI->obsImg1PositionSlider.RA;
    state->obsImg[1-1].DEC = formGUI->obsImg1PositionSlider.DEC;
    state->obsImg[2-1].RA  = formGUI->obsImg2PositionSlider.RA;
    state->obsImg[2-1].DEC = formGUI->obsImg2PositionSlider.DEC;
    state->obsImg[3-1].RA  = formGUI->obsImg3PositionSlider.RA;
    state->obsImg[3-1].DEC = formGUI->obsImg3PositionSlider.DEC;
    state->obsImg[4-1].RA  = formGUI->obsImg4PositionSlider.RA;
    state->obsImg[4-1].DEC = formGUI->obsImg4PositionSlider.DEC;
    state->obsImg[5-1].RA  = formGUI->obsImg5PositionSlider.RA;
    state->obsImg[5-1].DEC = formGUI->obsImg5PositionSlider.DEC;
    state->obsImg[6-1].RA  = formGUI->obsImg6PositionSlider.RA;
    state->obsImg[6-1].DEC = formGUI->obsImg6PositionSlider.DEC;
    state->obsImg[7-1].RA  = formGUI->obsImg7PositionSlider.RA;
    state->obsImg[7-1].DEC = formGUI->obsImg7PositionSlider.DEC;
    state->obsImg[8-1].RA  = formGUI->obsImg8PositionSlider.RA;
    state->obsImg[8-1].DEC = formGUI->obsImg8PositionSlider.DEC;

    state->obsImgIsVisible[1-1] = formGUI->obsImg1IsVisible;
    state->obsImgIsVisible[2-1] = formGUI->obsImg2IsVisible;
    state->obsImgIsVisible[3-1] = formGUI->obsImg3IsVisible;
    state->obsImgIsVisible[4-1] = formGUI->obsImg4IsVisible;
    state->obsImgIsVisible[5-1] = formGUI->obsImg5IsVisible;
    state->obsImgIsVisible[6-1] = formGUI->obsImg6IsVisible;
    state->obsImgIsVisible[7-1] = formGUI->obsImg7IsVisible;
    state->obsImgIsVisible[8-1] = formGUI->obsImg8IsVisible;

    DRAW_EXT_PARITY = formGUI->viewParity;
    DRAW_EXT_AS_POINTS = formGUI->viewExtSrcPoints;

    gridLimDistance = formGUI->gridLimDistance;
    SCALE_SPLANE = formGUI->tmpSrcScale/10;
    SCALE_IPLANE = formGUI->tmpImgScale/10;
}

void DrawSettingsForm(State *state, FormGUI *formGUI) {
    if (!formGUI->formActive) return;
    bool hasBeenUpdated = false;
    bool tmpCheckUpdate;

    Rectangle formRect = {panelUIposx + 10, panelUIposy + 110, panelUIsizx - 40, 30};

    int controlY = formRect.y + 20;
    int controlPadding = 35;
    float guiBoxSrcY = 0;

    if (dropdownLensModelEditMode || dropdownSourceModelEditMode) GuiLock();

    // LENS PARAMETERS --------------------------------------------------
    guiBoxSrcY = controlY;
    controlY += 10;
    Rectangle dropdownLensModelBounds = (Rectangle){formRect.x + 120, controlY, 90, 20};
    GuiLabel((Rectangle){formRect.x + 10, controlY, 100, 20}, "Lens Model:"); // Lens Mass
    // NOTE: GuiDropdownBox must draw after any other control that can be covered on unfolding
    // its at the end of this function
    controlY += controlPadding;

    if (formGUI->lensModelType != LTYPE_SPEMD) {
        // Lens mass parameter
        GuiLabel((Rectangle){formRect.x + 10, controlY, 100, 20}, "Mass :"); // Lens Mass (x10^12)
        tmpCheckUpdate = GuiSlider(
                            (Rectangle){formRect.x + 120, controlY, 90, 20},
                            TextFormat("%.2f", 0.1f),
                            TextFormat("%.2f", 2.0f),
                            &(formGUI->lensMassInput), 0.1f, 2.0f);
        hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
        controlY += controlPadding;
    }

    if (formGUI->lensModelType == LTYPE_SPEMD) {
        // Lens b parameter
        GuiLabel((Rectangle){formRect.x + 10, controlY, 100, 20}, "b (Ein rad):");
        tmpCheckUpdate = GuiSlider(
                            (Rectangle){formRect.x + 120, controlY, 90, 20},
                            TextFormat("%.2f", 0.1f),
                            TextFormat("%.2f", 3.0f),
                            &(formGUI->lensbInput), 0.1f, 3.0f);
        hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
        controlY += controlPadding;

        // Lens qh parameter
        GuiLabel((Rectangle){formRect.x + 10, controlY, 100, 20}, "qh (slope):");
        tmpCheckUpdate = GuiSlider(
                            (Rectangle){formRect.x + 120, controlY, 90, 20},
                            TextFormat("%.2f", 0.1f),
                            TextFormat("%.2f", 1.0f),
                            &(formGUI->lensqhInput), 0.1f, 1.0f);
        hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
        controlY += controlPadding;

        // // Lens rc parameter
        // GuiLabel((Rectangle){formRect.x + 10, controlY, 100, 20}, "rc :");
        // tmpCheckUpdate = GuiSlider(
        //                     (Rectangle){formRect.x + 120, controlY, 90, 20},
        //                     TextFormat("%.2f", 0.5f),
        //                     TextFormat("%.2f", 2.0f),
        //                     &(formGUI->lensrcInput), 1.0E-7, 1.0E-5);
        // hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
        // controlY += controlPadding;
    }

    if (formGUI->lensModelType != LTYPE_POINT) {
        // Lens f (axis ratio) parameter
        GuiLabel((Rectangle){formRect.x + 10, controlY, 100, 20}, "axis ratio f:");
        tmpCheckUpdate = GuiSlider(
                            (Rectangle){formRect.x + 120, controlY, 90, 20},
                            TextFormat("%.2f", 0.1f),
                            TextFormat("%.2f", 1.0f),
                            &(formGUI->lensFInput), 0.10f, 0.99f);
        hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
        controlY += controlPadding;
        // Lens PA (position angle) parameter
        GuiLabel((Rectangle){formRect.x + 10, controlY, 100, 20}, "PA [radians]:");
        tmpCheckUpdate = GuiSlider(
                            (Rectangle){formRect.x + 120, controlY, 90, 20},
                            TextFormat("%.2f", 0.0f),
                            TextFormat("%.2f", 3.14f),
                            &(formGUI->lensPAInput), 0, PI);
        hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
        controlY += controlPadding;
    }

    // Lens position controls
    GuiLabel((Rectangle){formRect.x + 10, controlY, 280, 20}, "Lens Position (as):");
    controlY += 20;
    GuiLabel((Rectangle){formRect.x + 55, controlY, 30, 20}, "RA:");
    tmpCheckUpdate = GuiSlider(
                        (Rectangle){formRect.x + 120, controlY, 90, 20},
                        TextFormat("%.2f", -1.0f),
                        TextFormat("%.2f",  1.0f),
                        &(formGUI->lensPositionSlider.RA), -1.0f, 1.0f);
    hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
    controlY += 20;
    GuiLabel((Rectangle){formRect.x + 55, controlY, 30, 20}, "DEC:");
    tmpCheckUpdate = GuiSlider(
                        (Rectangle){formRect.x + 120, controlY, 90, 20},
                        TextFormat("%.2f", -1.0f),
                        TextFormat("%.2f",  1.0f),
                        &(formGUI->lensPositionSlider.DEC), -1.0f, 1.0f);
    hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
    controlY += controlPadding;

    // TODO: implement external shear in all modes!
    GuiLabel((Rectangle){formRect.x + 10, controlY, 280, 20}, "Ext. Shear:");
    controlY += 20;
    GuiLabel((Rectangle){formRect.x + 55, controlY, 30, 20}, "g1:");
    tmpCheckUpdate = GuiSlider(
                        (Rectangle){formRect.x + 120, controlY, 90, 20},
                        TextFormat("%.1f", -0.50f),
                        TextFormat("%.1f", 0.50f),
                        &(formGUI->gamma_1), -0.50f, 0.50f);
    hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
    controlY += 20;
    GuiLabel((Rectangle){formRect.x + 55, controlY, 30, 20}, "g2:");
    tmpCheckUpdate = GuiSlider(
                        (Rectangle){formRect.x + 120, controlY, 90, 20},
                        TextFormat("%.1f", -0.50f),
                        TextFormat("%.1f", 0.50f),
                        &(formGUI->gamma_2), -0.50f, 0.50f);
    hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
    controlY += controlPadding;

    GuiGroupBox((Rectangle){formRect.x, guiBoxSrcY, formRect.width + 20, controlY - guiBoxSrcY - 10},
                            "Lens parameters");

    // SOURCE PARAMETERS --------------------------------------------------
    if (CURRENT_MODE == SRC_TO_LENS_MODE) {

        guiBoxSrcY = controlY;
        controlY += 10;

        // Source type (point or extended)
        GuiLabel((Rectangle){formRect.x + 10, controlY, 280, 20}, "Source Model:");
        Rectangle dropdownSourceModelBounds = (Rectangle){formRect.x + 120, controlY, 90, 20};
        // GuiCheckBox(
        //     (Rectangle){formRect.x + 120, controlY, 20, 20},
        //     "Extended",
        //     &(formGUI->useExtendedSource));
        controlY += controlPadding;

        // Source position controls
        GuiLabel((Rectangle){formRect.x + 10, controlY, 280, 20}, "Source Position (as):");
        controlY += 20;

        GuiLabel((Rectangle){formRect.x + 55, controlY, 30, 20}, "RA:");
        tmpCheckUpdate = GuiSlider(
                            (Rectangle){formRect.x + 120, controlY, 90, 20},
                            TextFormat("%.2f", -2.0f),
                            TextFormat("%.2f",  2.0f),
                            &(formGUI->sourcePositionSlider.RA), -2.0f, 2.0f);
        hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
        controlY += 20;
        GuiLabel((Rectangle){formRect.x + 55, controlY, 30, 20}, "DEC:");
        tmpCheckUpdate = GuiSlider(
                            (Rectangle){formRect.x + 120, controlY, 90, 20},
                            TextFormat("%.2f", -2.0f),
                            TextFormat("%.2f",  2.0f),
                            &(formGUI->sourcePositionSlider.DEC), -2.0f, 2.0f);
        hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
        controlY += controlPadding;

        // Source radius (only if extended)
        if (formGUI->sourceModelType == STYPE_EXT) {
            GuiLabel((Rectangle){formRect.x + 10, controlY, 100, 20}, "Radius (as):");
            tmpCheckUpdate = GuiSlider(
                            (Rectangle){formRect.x + 120, controlY, 90, 20},
                            TextFormat("%.2f", 0.05f),
                            TextFormat("%.2f", 1.0f),
                            &(formGUI->sourceRadius), 0.05f, 1.0f);
            hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
            controlY += controlPadding;

            tmpCheckUpdate = GuiCheckBox((Rectangle){formRect.x + 10, controlY, 20, 20} , "Points/Line", &(formGUI->viewExtSrcPoints));
            tmpCheckUpdate = GuiCheckBox((Rectangle){formRect.x + 120, controlY, 20, 20}, "Show parity", &(formGUI->viewParity));
            hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
            controlY += controlPadding;

        }

        GuiGroupBox((Rectangle){formRect.x, guiBoxSrcY, formRect.width + 20, controlY - guiBoxSrcY -10}, "Source parameters");

        // Dropbox for Source Model
        if (GuiDropdownBox(dropdownSourceModelBounds, dropdownSourceModelItems, &dropdownSourceModel, dropdownSourceModelEditMode))
            {
                dropdownSourceModelEditMode = !dropdownSourceModelEditMode;
                tmpCheckUpdate = true;
            }
        formGUI->sourceModelType = dropdownSourceModel;
        hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;

    }
    // OBSERVED IMAGE PARAMETERS --------------------------------------------------
    if (CURRENT_MODE == LENS_TO_SRC_MODE) {
        guiBoxSrcY = controlY;
        controlY += 10;
        // Source position controls
        GuiLabel((Rectangle){formRect.x + 70, controlY, 280, 20}, "Image Position (as):");
        controlY += 20;

        GuiLabel((Rectangle){formRect.x + 75, controlY, 30, 20}, "RA:");
        GuiLabel((Rectangle){formRect.x + 165, controlY, 30, 20}, "DEC:");
        controlY += 20;

        //for 1 to N_IMAGES observed images
        // TODO: Add color and visibility toggle for each observed image

        SkyCoord *imgPosTmp;
        bool *visObsImgTmp;
        for (int i = 0; i < 8; i++) {
            switch (i) {
                case 0:
                    imgPosTmp = &(formGUI->obsImg1PositionSlider);
                    visObsImgTmp = &(formGUI->obsImg1IsVisible);
                    break;
                case 1:
                    imgPosTmp = &(formGUI->obsImg2PositionSlider);
                    visObsImgTmp = &(formGUI->obsImg2IsVisible);
                    break;
                case 2:
                    imgPosTmp = &(formGUI->obsImg3PositionSlider);
                    visObsImgTmp = &(formGUI->obsImg3IsVisible);
                    break;
                case 3:
                    imgPosTmp = &(formGUI->obsImg4PositionSlider);
                    visObsImgTmp = &(formGUI->obsImg4IsVisible);
                    break;
                case 4:
                    imgPosTmp = &(formGUI->obsImg5PositionSlider);
                    visObsImgTmp = &(formGUI->obsImg5IsVisible);
                    break;
                case 5:
                    imgPosTmp = &(formGUI->obsImg6PositionSlider);
                    visObsImgTmp = &(formGUI->obsImg6IsVisible);
                    break;
                case 6:
                    imgPosTmp = &(formGUI->obsImg7PositionSlider);
                    visObsImgTmp = &(formGUI->obsImg7IsVisible);
                    break;
                case 7:
                    imgPosTmp = &(formGUI->obsImg8PositionSlider);
                    visObsImgTmp = &(formGUI->obsImg8IsVisible);
                    break;
                default:
                    imgPosTmp = &(formGUI->obsImg1PositionSlider);
                    visObsImgTmp = &(formGUI->obsImg1IsVisible);
                    break;
            }

            // DrawRectangleRec((Rectangle){formRect.x + 10, controlY, 20, 20}, colors[i % numColors]);
            DrawCircleV((Vector2){formRect.x + 20, controlY + 10}, 10, colors[i % numColors]);
            tmpCheckUpdate = GuiCheckBox((Rectangle){formRect.x + 38, controlY, 20, 20} , "", visObsImgTmp);
            hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
            tmpCheckUpdate = GuiSlider(
                                (Rectangle){formRect.x + 70, controlY, 80, 20},
                                "","",&(imgPosTmp->RA), -IMG_WIDTH_AS/2.0f, IMG_WIDTH_AS/2.0f);
            hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;

            tmpCheckUpdate = GuiSlider(
                                (Rectangle){formRect.x + 160, controlY, 80, 20},
                                "","",&(imgPosTmp->DEC), -IMG_HEIGHT_AS/2.0f, IMG_HEIGHT_AS/2.0f);
            hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
            controlY += 24;
        }
        controlY += 20;
        GuiGroupBox((Rectangle){formRect.x, guiBoxSrcY, formRect.width + 20, controlY - guiBoxSrcY - 10}, "Observed images");
    }

    // GuiDropdownBox closure --------------------------------------------------
    // NOTE: GuiDropdownBox must draw after any other control that can be covered on unfolding
    // Dropbox for Lens Model
    GuiUnlock();
    if (GuiDropdownBox(dropdownLensModelBounds, dropdownLensModelItems, &dropdownLensModel, dropdownLensModelEditMode))
        {
            dropdownLensModelEditMode = !dropdownLensModelEditMode;
            tmpCheckUpdate = true;
        }
    formGUI->lensModelType = dropdownLensModel;
    SELECTED_LENS_MODEL = dropdownLensModel;
    hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;

    if (hasBeenUpdated) { //Check if anything changed
        UpdateLensAndSourceModelsFromGUI(state, formGUI);
    }
}

typedef struct {
    float val;
    char *val_buffer;
    char *text_val;
    bool allow_edit;
} ValueBoxFloatControl;

ValueBoxFloatControl* CreateValueBoxFloatControl(float initialValue) {
    ValueBoxFloatControl *vbf = (ValueBoxFloatControl *)malloc(sizeof(ValueBoxFloatControl));
    vbf->val = initialValue;
    vbf->val_buffer = (char *)calloc(128, sizeof(char));
    vbf->text_val = (char *)calloc(128, sizeof(char));
    snprintf(vbf->text_val, 128, "%2.1f", vbf->val);
    vbf->allow_edit = false;
    return vbf;
}

void ImageFITSZoomGUIMenu(State *state, FormGUI *formGUI, GuiWindowFileDialogState* fileDialogState, ValueBoxFloatControl* vbf) {
    bool hasBeenUpdated = false;
    bool tmpCheckUpdate;
    bool editASValueBoxMode = false;

    Rectangle formRect = {ctrExtraPanX + 330, ctrExtraPanY - 120, 260, 30};
    float guiBoxSrcY = formRect.y;
    int controlY = formRect.y + 20;
    int controlPadding = 35;

    // Load Fits window --------------------------------------------------
    if (fileDialogState->windowActive) GuiLock();

    if (GuiButton((Rectangle){ formRect.x + 60, controlY, 140, 30}, GuiIconText(ICON_FILE_OPEN, "Open .fits"))){
        fileDialogState->windowActive = true;
    }
    GuiUnlock();
    controlY += controlPadding;

    // Set size of .fits in arcseconds -----------------------------------
    char textValue[RAYGUI_VALUEBOX_MAX_CHARS + 1] = "\0";
    float oldVal = vbf->val;
    if (GuiValueBoxFloat(
        (Rectangle){ formRect.x + 160, controlY, 40, 30},
        "side in as:",
        vbf->text_val,
        &vbf->val,
        vbf->allow_edit))
    {
        vbf->allow_edit = !vbf->allow_edit;
    }
    if (vbf->val < 0.1) { vbf->val = 0.1; }
    else if (vbf->val > 10.0) { vbf->val = 10.0; }
    snprintf(vbf->val_buffer, 128, "%2.1f", vbf->val);
    if (oldVal != vbf->val) {
        IMG_HEIGHT_AS = vbf->val;
        IMG_WIDTH_AS = vbf->val;
        UpdateImagePlane(state->imgPlane, img_fname, IMG_WIDTH_AS, IMG_HEIGHT_AS);
    }
    controlY += controlPadding;

    // SCALING PARAMETERS --------------------------------------------------
    controlY += 10;
    int paddingX = 15;
    GuiLabel((Rectangle){formRect.x + paddingX, controlY, 280, 20}, "Src plane");
    tmpCheckUpdate = GuiSlider(
                                (Rectangle){formRect.x + 120, controlY, 90, 20},
                                TextFormat("%.2f", 0.05f),
                                TextFormat("%.2f", 0.35f),
                                &(formGUI->tmpSrcScale), 0.05f, 0.35f);
    hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
    controlY += controlPadding;
    GuiLabel((Rectangle){formRect.x + paddingX, controlY, 280, 20}, "Img plane");
    tmpCheckUpdate = GuiSlider(
                                (Rectangle){formRect.x + 120, controlY, 90, 20},
                                TextFormat("%.2f", 0.05f),
                                TextFormat("%.2f", 0.35f),
                                &(formGUI->tmpImgScale), 0.05f, 0.35f);
    hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
    controlY += controlPadding;
    GuiCheckBox((Rectangle){formRect.x + paddingX, controlY, 20, 20},
                "Photometry",
                &(formGUI->viewFITS));
    GuiCheckBox((Rectangle){formRect.x + 120, controlY, 20, 20},
                formGUI->viewScaleBar ? "Scale" : "Scale",
                &(formGUI->viewScaleBar));
    controlY += controlPadding;

    if(CURRENT_MODE == SRC_TO_LENS_MODE) {
        GuiLabel((Rectangle){formRect.x + paddingX, controlY, 280, 20}, "Max grid dist.");
        tmpCheckUpdate = GuiSlider(
                                    (Rectangle){formRect.x + 120, controlY, 90, 20},
                                    TextFormat("%.2f", 0.0f),
                                    TextFormat("%.2f", 10.0f),
                                    &(formGUI->gridLimDistance), 0.0f, 10.0f);
        hasBeenUpdated = hasBeenUpdated || tmpCheckUpdate;
        controlY += controlPadding;
    }
    GuiGroupBox((Rectangle){formRect.x, guiBoxSrcY, formRect.width, controlY - formRect.y + 10}, "Load .fits & zoom");

    // GUI: Dialog Window
    GuiWindowFileDialog(fileDialogState);

    if (hasBeenUpdated) { //Check if anything changed
        UpdateLensAndSourceModelsFromGUI(state, formGUI);
    }
}


void UpdateFormGUI(State* state, FormGUI* formGUI){
    if (formGUI->formActive && state->lensMan->count > 0 && state->sourceMan->count > 0) {
        formGUI->lensModelType = state->lensMan->lenses[0].modelType;
        SELECTED_LENS_MODEL = formGUI->lensModelType;
        formGUI->lensMassInput = state->lensMan->lenses[0].M / 1e12;
        formGUI->lensbInput = state->lensMan->lenses[0].b;
        formGUI->lensqhInput = state->lensMan->lenses[0].qh;
        formGUI->lensrcInput = state->lensMan->lenses[0].rc;
        formGUI->lensFInput = state->lensMan->lenses[0].f;
        formGUI->lensPAInput = state->lensMan->lenses[0].PA;
        formGUI->lensPositionSlider.RA = state->lensMan->lenses[0].coord.RA;
        formGUI->lensPositionSlider.DEC = state->lensMan->lenses[0].coord.DEC;
        formGUI->gamma_1 = state->gamma1;
        formGUI->gamma_2 = state->gamma2;

        formGUI->sourceModelType = state->sourceMan->sources[0].modelType;
        formGUI->sourcePositionSlider.RA = state->sourceMan->sources[0].coord.RA;
        formGUI->sourcePositionSlider.DEC = state->sourceMan->sources[0].coord.DEC;

        formGUI->obsImg1PositionSlider.RA  = state->obsImg[1-1].RA;
        formGUI->obsImg1PositionSlider.DEC = state->obsImg[1-1].DEC;
        formGUI->obsImg2PositionSlider.RA  = state->obsImg[2-1].RA;
        formGUI->obsImg2PositionSlider.DEC = state->obsImg[2-1].DEC;
        formGUI->obsImg3PositionSlider.RA  = state->obsImg[3-1].RA;
        formGUI->obsImg3PositionSlider.DEC = state->obsImg[3-1].DEC;
        formGUI->obsImg4PositionSlider.RA  = state->obsImg[4-1].RA;
        formGUI->obsImg4PositionSlider.DEC = state->obsImg[4-1].DEC;
        formGUI->obsImg5PositionSlider.RA  = state->obsImg[5-1].RA;
        formGUI->obsImg5PositionSlider.DEC = state->obsImg[5-1].DEC;
        formGUI->obsImg6PositionSlider.RA  = state->obsImg[6-1].RA;
        formGUI->obsImg6PositionSlider.DEC = state->obsImg[6-1].DEC;
        formGUI->obsImg7PositionSlider.RA  = state->obsImg[7-1].RA;
        formGUI->obsImg7PositionSlider.DEC = state->obsImg[7-1].DEC;
        formGUI->obsImg8PositionSlider.RA  = state->obsImg[8-1].RA;
        formGUI->obsImg8PositionSlider.DEC = state->obsImg[8-1].DEC;

        formGUI->obsImg1IsVisible = state->obsImgIsVisible[1-1];
        formGUI->obsImg2IsVisible = state->obsImgIsVisible[2-1];
        formGUI->obsImg3IsVisible = state->obsImgIsVisible[3-1];
        formGUI->obsImg4IsVisible = state->obsImgIsVisible[4-1];
        formGUI->obsImg5IsVisible = state->obsImgIsVisible[5-1];
        formGUI->obsImg6IsVisible = state->obsImgIsVisible[6-1];
        formGUI->obsImg7IsVisible = state->obsImgIsVisible[7-1];
        formGUI->obsImg8IsVisible = state->obsImgIsVisible[8-1];

        formGUI->sourceRadius = state->sourceMan->sources[0].radius;

        dropdownLensModel = formGUI->lensModelType;
        dropdownSourceModel = formGUI->sourceModelType;
    }
}

// Add this to initialize the form state
void InitFormGUI(State* state, FormGUI* formGUI) {
    if (state->lensMan->count > 0 && state->sourceMan->count > 0) {
        formGUI->formActive = true;

        formGUI->lensModelType = state->lensMan->lenses[0].modelType;
        formGUI->lensMassInput = state->lensMan->lenses[0].M / 1e12;
        formGUI->lensbInput = state->lensMan->lenses[0].b;
        formGUI->lensqhInput = state->lensMan->lenses[0].qh;
        formGUI->lensrcInput = state->lensMan->lenses[0].rc;
        formGUI->lensFInput = state->lensMan->lenses[0].f;
        formGUI->lensPAInput = state->lensMan->lenses[0].PA;
        formGUI->lensPositionSlider.RA = state->lensMan->lenses[0].coord.RA;
        formGUI->lensPositionSlider.DEC = state->lensMan->lenses[0].coord.DEC;
        formGUI->gamma_1 = state->gamma1;
        formGUI->gamma_2 = state->gamma2;

        formGUI->sourceModelType = state->sourceMan->sources[0].modelType;
        formGUI->sourcePositionSlider.RA = state->sourceMan->sources[0].coord.RA;
        formGUI->sourcePositionSlider.DEC = state->sourceMan->sources[0].coord.DEC;

        formGUI->obsImg1PositionSlider.RA  = state->obsImg[1-1].RA;
        formGUI->obsImg1PositionSlider.DEC = state->obsImg[1-1].DEC;
        formGUI->obsImg1PositionSlider.RA  = state->obsImg[2-1].RA;
        formGUI->obsImg1PositionSlider.DEC = state->obsImg[2-1].DEC;
        formGUI->obsImg1PositionSlider.RA  = state->obsImg[3-1].RA;
        formGUI->obsImg1PositionSlider.DEC = state->obsImg[3-1].DEC;
        formGUI->obsImg1PositionSlider.RA  = state->obsImg[4-1].RA;
        formGUI->obsImg1PositionSlider.DEC = state->obsImg[4-1].DEC;
        formGUI->obsImg1PositionSlider.RA  = state->obsImg[5-1].RA;
        formGUI->obsImg1PositionSlider.DEC = state->obsImg[5-1].DEC;
        formGUI->obsImg1PositionSlider.RA  = state->obsImg[6-1].RA;
        formGUI->obsImg1PositionSlider.DEC = state->obsImg[6-1].DEC;
        formGUI->obsImg1PositionSlider.RA  = state->obsImg[7-1].RA;
        formGUI->obsImg1PositionSlider.DEC = state->obsImg[7-1].DEC;
        formGUI->obsImg1PositionSlider.RA  = state->obsImg[8-1].RA;
        formGUI->obsImg1PositionSlider.DEC = state->obsImg[8-1].DEC;

        formGUI->obsImg1IsVisible = state->obsImgIsVisible[1-1];
        formGUI->obsImg2IsVisible = state->obsImgIsVisible[2-1];
        formGUI->obsImg3IsVisible = state->obsImgIsVisible[3-1];
        formGUI->obsImg4IsVisible = state->obsImgIsVisible[4-1];
        formGUI->obsImg5IsVisible = state->obsImgIsVisible[5-1];
        formGUI->obsImg6IsVisible = state->obsImgIsVisible[6-1];
        formGUI->obsImg7IsVisible = state->obsImgIsVisible[7-1];
        formGUI->obsImg8IsVisible = state->obsImgIsVisible[8-1];

        formGUI->sourceRadius = state->sourceMan->sources[0].radius;

        dropdownLensModel = formGUI->lensModelType;
        dropdownSourceModel = formGUI->sourceModelType;

    } else {
        formGUI->formActive = true;

        dropdownLensModel = LTYPE_SPEMD;
        formGUI->lensModelType = LTYPE_SPEMD;
        formGUI->lensMassInput = LENS_MASS/1e12;
        formGUI->lensbInput = 1.0f;
        formGUI->lensqhInput = 0.5f;
        formGUI->lensrcInput = 1.0000000000E-6f;
        formGUI->lensFInput = 0.3f;
        formGUI->lensPAInput = 0.f;
        formGUI->lensPositionSlider   = (SkyCoord){0.0f, 0.0f};
        formGUI->gamma_1   = 0.0f;
        formGUI->gamma_2   = 0.0f;

        dropdownSourceModel = STYPE_EXT;
        formGUI->sourceModelType = STYPE_EXT;
        formGUI->sourceRadius = 0.5f;
        formGUI->sourcePositionSlider = (SkyCoord){0.5f, 0.5f};

        formGUI->obsImg1PositionSlider  = (SkyCoord){0.0f, 0.0f};
        formGUI->obsImg2PositionSlider  = (SkyCoord){0.0f, 0.0f};
        formGUI->obsImg3PositionSlider  = (SkyCoord){0.0f, 0.0f};
        formGUI->obsImg4PositionSlider  = (SkyCoord){0.0f, 0.0f};
        formGUI->obsImg5PositionSlider  = (SkyCoord){0.0f, 0.0f};
        formGUI->obsImg6PositionSlider  = (SkyCoord){0.0f, 0.0f};
        formGUI->obsImg7PositionSlider  = (SkyCoord){0.0f, 0.0f};
        formGUI->obsImg8PositionSlider  = (SkyCoord){0.0f, 0.0f};

        formGUI->obsImg1IsVisible = 0;
        formGUI->obsImg2IsVisible = 0;
        formGUI->obsImg3IsVisible = 0;
        formGUI->obsImg4IsVisible = 0;
        formGUI->obsImg5IsVisible = 0;
        formGUI->obsImg6IsVisible = 0;
        formGUI->obsImg7IsVisible = 0;
        formGUI->obsImg8IsVisible = 0;
    }

    formGUI->gridLimDistance = gridLimDistance;
    formGUI->tmpSrcScale = SCALE_SPLANE * 10;
    formGUI->tmpImgScale = SCALE_IPLANE * 10;
    formGUI->viewFITS = true;
    formGUI->viewScaleBar = true;
    formGUI->viewParity = DRAW_EXT_PARITY;
    formGUI->viewExtSrcPoints = DRAW_EXT_AS_POINTS;

    TraceLog(LOG_DEBUG, "FormGUI initialized");
}

void FreeFormGUI(FormGUI* formGUI) {
    if (formGUI != NULL) {
        free(formGUI);
    }
}

void MenuImGUI(State* state, FormGUI* formGUI) {
    GuiPanel((Rectangle){panelUIposx, panelUIposy, panelUIsizx, panelUIsizy}, "Menu");
    if (GuiButton((Rectangle){panelUIposx + 20, panelUIposy + 40, panelUIsizx - 40, 30},
        Get_Mode_Text(CURRENT_MODE))) {
            CURRENT_MODE = (CURRENT_MODE + 1) % NUM_MODES;
    }

    // UI elements to make new random image
    if (GuiButton((Rectangle){panelUIposx + 20, panelUIposy + 80, panelUIsizx - 40, 30},
        "Generate random lens - source config.")) {
        // Generate a random position for the source in the source plane
        ReplaceLastLensWithRandomLens(state->lensMan);
        ReplaceLastSourceWithRandomSource(state->sourceMan);
        state->gamma1 = 0.0;
        state->gamma2 = 0.0;
        UpdateLensAndSourceModels(state);
        UpdateFormGUI(state, formGUI);
        TraceLog(LOG_DEBUG, "Updated state after generating random lens-source");
    }

    // if (GuiButton((Rectangle){panelUIposx + 20, panelUIposy + 110, panelUIsizx - 40, 30},
    //     "Edit Lens & Source")) {
    //     formGUI->formActive = !formGUI->formActive;
    //     UpdateFormGUI(state, formGUI);
    // }

    // UI elements to toggle DARKMODE
    if (GuiButton((Rectangle){panelUIposx + 20, panelUIposy + panelUIsizy - 40, panelUIsizx - 40, 30},
        DARKMODE ?
        "Darkmode: ON " :
        "Darkmode: OFF")) {
            DARKMODE = !DARKMODE;
            DARKMODE ? GuiLoadStyleAmber() : GuiLoadStyleDefault();
    }
}

//------------------------------------------------------------------------------------
// Program main entry point
//------------------------------------------------------------------------------------
int main(void)
{
    // Initialization
    //--------------------------------------------------------------------------------------
    InitWindow(screenWidth, screenHeight, "lensfocus - quickly explore lens models");
    SetTraceLogLevel(CURRENT_LOG_MODE);

    DARKMODE ? GuiLoadStyleAmber() : GuiLoadStyleDefault();

    //------------------------------------------------------------------------------------
    // Custom file dialog
    GuiWindowFileDialogState fileDialogState = InitGuiWindowFileDialog(GetWorkingDirectory());
    bool exitWindow = false;
    strcpy(img_fname, default_img_fname);

    // Init default ValueBoxFloatControl for .fits size
    ValueBoxFloatControl* vbf = CreateValueBoxFloatControl((float)IMG_WIDTH_AS);

    // Init state
    CURRENT_MODE = LENS_TO_SRC_MODE;

    State* state  = StateInit();
    int lens_id   = AddLens(state->lensMan, GetRandomLens());
    int source_id = AddSource(state->sourceMan, GetRandomSource());

    UpdateLensAndSourceModels(state);

    // Init form state
    FormGUI* formGUI = (FormGUI*)malloc(sizeof(FormGUI));
    InitFormGUI(state,formGUI);

    SetTargetFPS(60);                // Set our app to run at 60 frames-per-second
    //--------------------------------------------------------------------------------------

    // Main loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        if (fileDialogState.SelectFilePressed)
        {
            // Load image file (if supported extension)
            if (IsFileExtension(fileDialogState.fileNameText, ".fits"))
            {
                strcpy(img_fname, TextFormat("%s" PATH_SEPERATOR "%s", fileDialogState.dirPathText, fileDialogState.fileNameText));
                UpdateImagePlane(state->imgPlane, img_fname, IMG_WIDTH_AS, IMG_HEIGHT_AS);
            }

            fileDialogState.SelectFilePressed = false;
        }

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();
        SetupSrcImgPlanes();
        MenuImGUI(state, formGUI);
        DrawSettingsForm(state, formGUI);
        //----------------------------------------------------------------------------------
        if(formGUI->viewFITS) {
            double scale = (IMG_WIDTH_AS / SCALE_IPLANE);
            int width_pix  = state->imgPlane->textureImageData->width;
            int height_pix = state->imgPlane->textureImageData->height;
            Rectangle srcRect = {
                0, 0,
                width_pix,
                height_pix
            };
            Rectangle dstRect = {
                ctrImgPlaneX - scale/2,
                ctrImgPlaneY - scale/2,
                scale,
                scale
            };
            Vector2 origin = {0, 0};
            float rotation = 0.0f;
            DrawTexturePro(*state->imgPlane->textureImageData,
                            srcRect, dstRect, origin, rotation, WHITE);
        }
        switch (CURRENT_MODE) {
            case LENS_TO_SRC_MODE: // Mode: Trace rays from lens to source plane
                if (state->lensMan->count>0 && state->sourceMan->count>0){
                    if (formGUI->viewScaleBar) {
                        DrawScaleBar();
                    }
                    DrawAllObsImages(state);
                    DrawCenterPaLens(state);
                    DrawAllCausticAndCriticalLine(state);
                }
                break;
            case SRC_TO_LENS_MODE: // Mode: Simple Lens Equation solver
                if (state->lensMan->count>0 && state->sourceMan->count>0){
                    if (formGUI->viewScaleBar) {
                        DrawScaleBar();
                    }
                    DrawAllSources(state);
                    DrawAllImages(state);
                    DrawAllCausticAndCriticalLine(state);
                }
                break;
            case NUM_MODES: // Mode: NUM_MODES
                // This mode is only to get the size of the enum
                // Do nothing here.
                break;
            default: // Mode: Unknown
                FreeFormGUI(formGUI);
                StateDestroy(state);
                CloseWindow();
                return -1; // TODO: Implement this mode;
        }

        //----------------------------------------------------------------------------------
        DisplayStateInfo(state);
        ImageFITSZoomGUIMenu(state, formGUI, &fileDialogState, vbf);
        EndDrawing();
    }

    FreeFormGUI(formGUI);
    StateDestroy(state);
    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
