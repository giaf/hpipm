SOURCES = [ 'qp_data.c ', ...
            'hpipm_sfun.c ', ...
          ];

HPIPM_MAIN_FOLDER = getenv('HPIPM_MAIN_FOLDER');
BLASFEO_MAIN_FOLDER = getenv('BLASFEO_MAIN_FOLDER');

INCS = [ ' -I', HPIPM_MAIN_FOLDER, '/include/ ', ...
         ' -I', BLASFEO_MAIN_FOLDER, '/include/ '];

CFLAGS  = ' -O';

LIB_PATH = [ ' -L', HPIPM_MAIN_FOLDER, '/lib/ ', ...
         ' -L', BLASFEO_MAIN_FOLDER, '/lib/ '];

LIBS = '-lhpipm -lblasfeo -lm'; 
    
eval( [ 'mex -v -output  hpipm_solver ', ...
    CFLAGS, INCS, ' ', SOURCES, LIB_PATH, ' ', LIBS ]);

disp( [ 'hpipm_solver', '.', ...
    eval('mexext'), ' successfully created!'] );

