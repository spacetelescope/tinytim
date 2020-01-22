#define  VERSION_NUM  7.5

/* 1-8 = WF/PC-1 */
#define  FOC_F48      9
#define  FOC_F96     10
#define  WFPC2_PC    11
#define  WFPC2_WFC_2 12
#define  WFPC2_WFC_3 13
#define  WFPC2_WFC_4 14
#define  FOC2_F48    15
#define  FOC2_F96    16
#define  NICMOS_1    17
#define  NICMOS_2    18
#define  NICMOS_3    19
#define  STIS_CCD    20
#define  STIS_NUV    21
#define  STIS_FUV    22
#define  ACS_WFC1    23
#define  ACS_WFC2    24
#define  ACS_HRC     25
#define  ACS_HRC_OFFSPOT     26
#define  ACS_SBC     27
#define  WFC3_UVIS1    28
#define  WFC3_UVIS2    29
#define  WFC3_IR     30
#define  NEWNICMOS_1 31
#define  NEWNICMOS_2 32
#define  NEWNICMOS_3 33

#define  M_PI  3.14159265358979323846
#define  DIAM_HST_MICRONS    2.4e6
#define  MAX_WAVES      200
#define  MAX_SPECTRA     50
#define  MAX_POSITIONS  1000
#define  LASTZERNIKE    22
#define  MAX_INSTR      30
#define  NUM_ACS_KERNELS 3
#define  NUM_WFC3_KERNELS 3
#define  MAX_STRING  550
#define  MAX_PSFS       1000

#define  MAX_THREADS 5

#define  LIMITSIZE   49

#define  KX  3
#define  KY  3
#define  KZ4 2

int Num_threads;

enum filetypes { MAP_FILE, PSF_FILE, DISTORTED_PSF_FILE, IMAGE_FILE };

struct paramstruct {
    float   ota_secondary_radius,
        ota_spider_width,
        ota_pad_v2[3],
        ota_pad_v3[3],
        ota_pad_radius[3],

        camera_entrance_radius,
        camera_secondary_radius,
        camera_spider_width,

        min_wavelength, max_wavelength,
        pixel_size,	/* normal detector pixel size in arcsec */
        psf_scale, 	/* pixel size of PSF produced by tiny2 in arcsec */
        psf_size,	/* diameter of PSF produced by tiny2 in arcsec */

        v2, v3,   	/* OTA axial location of center of detector */
        v2c, v3c,       /* location of camera (not detector) center */
        x_ref, y_ref,   /* reference location for ACS transforms */
 
        /* transformation coefficients for undistorted X,Y to detector X,Y */

        xcyc_to_x[14], xcyc_to_y[14],

        /* transformation coefficients for detector X,Y to undistorted X,Y */

        xy_to_xc[14], xy_to_yc[14],

        theta,
        caber[5][5],

        wfc_vig_radius, wfc_vig_xc, wfc_vig_yc, wfc_vig_inside,
        obsc_center_x, obsc_center_y,

        kernel[3][3],
        weighted_kernel[3][3],

                nic_pad_v2[3],
                nic_pad_v3[3],
                nic_pad_x[3],
                nic_pad_y[3],
                nic_pad_rot[3],
        pam_to_z4,
        z5_vs_z4[KZ4],z6_vs_z4[KZ4],z7_vs_z4[KZ4],z8_vs_z4[KZ4],z11_vs_z4[KZ4],
        z4_vs_xy[KY][KX], z5_vs_xy[KY][KX], z6_vs_xy[KY][KX], z7_vs_xy[KY][KX],
        z8_vs_xy[KY][KX],

        acs_outer_radius, acs_secondary_radius, acs_spider_width,
        acs_pad_radius, acs_stop_shift, acs_stop_angle,
        acs_kernel_wavelength[NUM_ACS_KERNELS],
        acs_kernel[NUM_ACS_KERNELS][3][3],
        kernel_xy[3][3],
        acs_aber[5][6][6], blur_xy_sigma[10][6][6], blur_xy_wavelength[10],
        undistorted_scale,

        stis_stop_x_offset, stis_stop_y_offset,
        stis_stop_x_radius, stis_stop_y_radius,

        wfc3_mask_x, wfc3_mask_y,
                wfc3_pad_v2[3],
                wfc3_pad_v3[3],
                wfc3_pad_x[3],
                wfc3_pad_y[3],
                wfc3_pad_rot[3],
                wfc3_aber[5][6][6],
        wfc3_kernel_wavelength[NUM_WFC3_KERNELS],
        wfc3_kernel[NUM_WFC3_KERNELS][3][3];

    char	tt3_param_file[MAX_STRING], scene_input_file[MAX_STRING],
        scene_output_file[MAX_STRING], scene_psf_file[MAX_STRING];
    float	scene_pixel_scale, scene_psf_x[MAX_PSFS], scene_psf_y[MAX_PSFS],
        scene_psf_flux[MAX_PSFS];
    int	scene_x_center, scene_y_center, num_scene_psfs;
    int	blur_xy_num_wavelengths;

    float   crit_scale[MAX_WAVES],
        wavelength[MAX_WAVES], weight[MAX_WAVES], wave_mag,
        jitter_major, jitter_minor, jitter_angle,
        sampling, zval[30], fnumber,
        best_focus_pam, observation_pam, ebmv, av,
        focus, xastig, yastig, xcoma, ycoma, spherical,
        xclover, yclover,
        nicmos_mask_x_off, nicmos_mask_y_off, nicmos_mask_rotation,
        ota_offset, weight_limit;
    int     chip, crit_dim[MAX_WAVES], current_pos, int_dim, num_waves,
        num_pos, use_map, write_crit_psf, write_pupil, write_wave,
        zset[30], scatter_flag, distorted_dim, sub_factor,
        adjust_for_focus, adjust_for_xy, smart_skip, add_halo,
        x[MAX_POSITIONS], y[MAX_POSITIONS], num_comments, old_nicmos;
    char    pupil_name[MAX_STRING], zernike_name[MAX_STRING],
        filter_name[MAX_STRING], root_name[MAX_STRING], camera_name[50],
        filter_file[MAX_STRING], spectrum_file[MAX_STRING],
        *comment[40];
} Pars;

struct DateStruct {
    int	day, month, year;
};

typedef struct Complex { float r, i; } complex;

#define  FLOATING  1
#define  COMPLEX   2

#define  ERROR     0
#define  OKAY      1

#define  TRUE      1
#define  FALSE     0

#define  RIGHT_JUSTIFY  1
#define  LEFT_JUSTIFY   0
#define  NO_VALUE      -1

#define  FITS_FILE      0
#define  STSDAS_FILE    1

/* acs.c */
float Compute_kernel_xy( int xdet, int ydet, float lambda );
void Compute_acs_kernel( float wavelength, float weight );
void Compute_wfc3_kernel( float wavelength, float weight );
void Convolve_kernel( float kernel1[3][3], float kernel2[3][3] );
float **Read_psf( int position );
float **Read_scene( int *nx, int *ny );
float **Convolve_scene( float **scene, int *scene_nx, int *scene_ny, float **psf_in );
complex **Resample( float **image, int nx, int ny, int dim_out, float mag );
complex **Copy_to_center( float **image, int nx, int ny, int dim_out );
void Read_optional_parameters( void );
void Add_psfs( float **psf, float **image, int nx, int ny );

/* compgrid.c */
void Compute_grid_dims( float *psf_diameter, int chip, 
    int num_waves, float *waves, int *crit_dim, int *int_dim );

/* convotf.c */
void Convolve_otf( complex **crit_psf, int dim, float data_scale, 
    float pixel_scale );

/* cpupil.c */
void Compute_pupil( int dim, int psf_x, int psf_y, float **pupil );

/* distort.c */
void Init_damped_interp( void );
/* float Damped_interp( float **image, float x, float y, int nx, int ny, int thread ); */
float Damped_interp( float **image, float x, float y, int nx, int ny );
void Distort( float **in, int nx, int ny, int x_det, int y_det, float **out, int out_nx, int out_ny, int camera );

/* entry.c */
void Init_entry_list( void );
void Delete_entry_list( void );
void Write_entry_list( FILE *file );
void Store_entry( char *entry );
char *Get_string( FILE *file );
char *Get_entry( FILE *file );
void Hack_z4_entry(float waves_secondary_mirror_despace);

/* fft.c */
void fft2d( complex *image, int dim, int direction );

/* fitsio.c */
float **Read_FITS( char *filename, int *nx, int *ny );
void Write_FITS( char *filename, float *image, int nx, int ny, 
    int type, int overwrite );

/* halo.c */
void Halo_f1042m( float **image, int dim, float subsamp );
void Hrc_halo( float **psf, float wavelength, int psf_x, int psf_y );
float Acs_pixel_area( int x, int y, int camera );

/* image.c */
void Report_num_alloc( void );
float **Alloc_image( int nx, int ny );
void Free_image( float **image );
complex **Alloc_complex_image( int nx, int ny );
void Free_complex_image( complex **image );
void Shift_to_center( float **image, int dim );
void Shift_to_origin( float **image, int dim );
void Shift_complex_to_center( complex **image, int dim );
void Shift_complex_to_origin( complex **image, int dim );

/* inputs.c */
char *Find_string( char *string1a, char *string2a );
void Str_up_case( char *string );
char *Str_token( char *string, char *token );
void Str_get_floats( char *string, int num_floats, float *floats );
int Get_filter( void );
int Read_positions( char *filename, int *x, int *y );
int Read_table( char *filename, float *wave, float *weight );

/* inter.c */
void Interpolate( float **old_image, int old_dim, float old_start_x, 
    float old_start_y, float old_spacing, float **new_image,
    int new_dim, float new_start_x,
    float new_start_y, float new_spacing );

/* intpsf.c */
void Integrate_psf( complex **psf_in, float scale_in, int size_in,
        float **psf_out, float scale_out, int size_out, 
    float jitter_major, float jitter_minor, float jitter_angle );

/* jitter.c */
void Compute_jitter( complex **crit_psf, int dim, 
    float wavelength, float jitter );

/* map.c */
float **Compute_map( int map_dim, int chip, float theta, 
    float v2, float v3 );

/* misc.c */
int roundval( float value );
void Convolve( float **image, float kernel[3][3], int nx, int ny );
void Flux_normalize( float *array, int nx, int ny );
float sinc( float x );
void Str_trim( char *string );
int Swap_test( void );
void strcompress( char *string );
void Get_keyword_value( char *entry, char *keyword, char *value );
void Swap_int( char *data, int number );
void Min_max( float *data, int num_pix, float *min_val, float *max_val );
void xstrncpy( char *dest, char *source, int number );
int Write_image( char *filename, void *image, int nx, int ny, int data_type );
int Read_image( char *filename, void *image, int nx, int ny, int data_type );

/* monopsf.c */
void Compute_mono_psf( int w, float **mono_psf, complex **pupil );

/* opd.c */
void Compute_opd( int dim, float **opd, int psf_x, int psf_y );

/* outputs.c */
void Write_params( char *param_file, struct DateStruct *obs_date );

/* param.c */
void Read_params( char *param_file );

/* polypsf.c */
void Compute_poly_psf( int pos_index, float **poly_psf, float **mono_psf );

/* psf.c */
void Compute_crit_psf( complex **pupil, int dim );

/* pupil.c */
void Draw_ellipse( float **image, int nx, int ny, float xc, float yc,
    float x_radius, float y_radius );
void Smooth( float **image, int nx, int ny, int box_size );
void Dark_inside_circle( float **image, int nx, int ny, 
    float xc, float yc, float radius );
void Dark_outside_circle( float **image, int nx, int ny, 
    float xc, float yc, float radius );
void Clear_inside_circle( float **image, int nx, int ny, 
    float xc, float yc, float radius );
void Dark_rectangle( float **image, int nx, int ny, 
    float xc, float yc, float x_length, float y_length );
void Dark_rot_rect( float **image, int nx, int ny, float xc, float yc,
        float x_length, float y_length, float theta );

/* rdpupil.c */
void Read_pupil_data( int camera, FILE *file );
void Read_pupil_table( int camera, char *table_name );
void Read_geom_data( FILE *file );

/* rotate.c */
float **Rotate_image( float **image, int nx, int ny, 
    float xc, float yc, float theta );
float **Rotate( float **image, int old_nx, int old_ny, float xc, float yc,
                int new_nx, int new_ny, float new_xc, float new_yc,
                float theta, float mag );

/* sincint.c */
void Sinc_int2d( float **old_image, int old_dim, float **new_image, 
    int new_dim, float old_start, float old_scale, float new_start,
    float new_scale, float convergence, float delta, float pixel_size,
        float limit_radius, int limit_size );

/* spectrum.c */
void Get_spectrum( void );

/* system.c */
void Delete_file( char *filename );
void Default_dir( char *string, char *filename );

/* v2v3.c */
void XY_to_V2V3( float x, float y, float *v2, float *v3, int camera );
void V2V3_to_XY( float *x, float *y, float v2, float v3, int camera );
void v2v3( int camera, int x, int y, float *v2, float *v3 );

/* vignet.c */
void Wfc_vignette( float **image, int nx, int ny );

/* zernike.c */
void Read_zernike_data( FILE *file );
void Get_aberrations( struct DateStruct *obs_date );
void Compute_aberrations( int psf_x, int psf_y );

