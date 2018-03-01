//---------------------------------------------------------------------------

#ifndef ufftH
#define ufftH
//---------------------------------------------------------------------------


void complex_16_express ( DWORD v,
                          short &r,
                          short &i );

DWORD complex_16_abs ( DWORD v );

void fft (   DWORD  T[], DWORD  F[],
             DWORD II[], DWORD IO[],
           LPDWORD  C[], size_t N );

void init_fft_coff ( LPDWORD fc[], size_t N );
void gen_w ( DWORD W[], size_t N );

LPDWORD complex_16_array_create ( size_t N );

void complex_16_array_delete ( LPDWORD );
void fft_shuf_t ( DWORD ti[], DWORD to[], size_t n );
void fft_shuf_f ( DWORD fi[], DWORD fo[], size_t N );

LPDWORD * fft_coff_create ( size_t N );
     void fft_coff_delete ( LPDWORD * coff, size_t N );

size_t integer_log2 ( size_t n );


struct s_fcomplex;
typedef struct s_fcomplex fcomplex_t, *fcomplex_p;

struct s_fcomplex
{
  float m_r;
  float m_i;
  static void add ( fcomplex_p c1,
                    fcomplex_p c2,
                    fcomplex_p r );
  static void sub ( fcomplex_p c1,
                    fcomplex_p c2,
                    fcomplex_p r );
  static void mul ( fcomplex_p c1,
                    fcomplex_p c2,
                    fcomplex_p r );
  static void mul2 ( fcomplex_p c1,
                    fcomplex_p c2,
                    fcomplex_p r );
  static void mul_add( fcomplex_p s, fcomplex_p c,
                       fcomplex_p r );
  static void mul_sub( fcomplex_p s, fcomplex_p c,
                       fcomplex_p r );
  static void mul_add_sub( fcomplex_p s,  fcomplex_p c,
                           fcomplex_p ra, fcomplex_p rs );
  static void mul_add_sub2( fcomplex_p s,  fcomplex_p c,
                           fcomplex_p psum, fcomplex_p pdef );
  static void load_from_IQ16 ( fcomplex_p r, DWORD IQ[], size_t N );
  static void load_from_IQ16_ex ( fcomplex_p r, DWORD IQ[], size_t N );
  static fcomplex_p array_create ( size_t N );
  static void array_delete ( fcomplex_p array );
  void set  ( float r, float i );
  void set  ( DWORD IQ );
  s_fcomplex ( float r, float i );
  s_fcomplex ( float r );
  s_fcomplex ( void );
};

class TFFT
{
private:
  size_t m_N;
  fcomplex_p *m_c;
  fcomplex_p m_ii;
  fcomplex_p m_io;
protected:
  static size_t fix_N_if_error ( size_t a_N );
  void set ( size_t N );
  void clr ( void );
  static void coff_init ( fcomplex_p c[], size_t N );
  static fcomplex_p * coff_create ( size_t N );
  static void coff_delete ( fcomplex_p c[], size_t N );
  void round ( fcomplex_p S, fcomplex_p C, fcomplex_p R );
  void round3 ( fcomplex_p S, fcomplex_p C, fcomplex_p R );
  void round2 ( fcomplex_p S, fcomplex_p C, fcomplex_p R );
  void shuf_t ( fcomplex_p ti, fcomplex_p to );
  void shuf_f ( fcomplex_p fi, fcomplex_p fo );
public:
  TFFT ( size_t a_N );
  ~TFFT ();
  void fft ( fcomplex_p t, fcomplex_p f );
};
#endif
