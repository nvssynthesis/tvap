#include "m_pd.h"
#include <math.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif
#define PI 3.141592653589793

/* ------------------------ tvap~ ----------------------------- */

/* time varying allpass 
2 arguments:
    -f_pi, the frequency of pi phase shift
    -f_b,  the bandwidth
2 inlets:
    -x[n], the signal to filter
        -also accepts f_b as message
    -f_pi[n], which also can be a signal
        -also accepts f_pi as message

TODO: 
	-create tangent and cosine table instead of tan(), sin(), cos()
	-maybe it could be useful to convert bandwidth to Q
*/
typedef struct tvap_ctl 
{
    t_sample z1, z2;
    t_float fs;
} t_tvap_ctl;

typedef struct _tvap_tilde
{
    t_object x_obj; 
    t_float x_sr; 
    t_float x_f_pi;
    t_float x_f_b; 
    t_tvap_ctl x_cspace;
    t_tvap_ctl *x_ctl;
    t_float x_f;
} t_tvap_tilde;

t_class *tvap_tilde_class;

static void tvap_set_f_pi(t_tvap_tilde *x, t_floatarg f_pi);

static void *tvap_tilde_new(t_floatarg f_pi, t_floatarg f_b)
{
    t_tvap_tilde *x = (t_tvap_tilde *)pd_new(tvap_tilde_class);
    signalinlet_new(&x->x_obj, f_pi);
    signalinlet_new(&x->x_obj, f_b);
    outlet_new(&x->x_obj, &s_signal);

    x->x_sr = 44100; // just make it something, then dsp corrects it
    x->x_ctl = &x->x_cspace; // effectively initializes x_ctl to 0's
    x->x_cspace.z1 = 0;
    x->x_cspace.z2 = 0;
    tvap_set_f_pi(x, f_pi);
    x->x_f = 0; 

    return (x);
}

float f_pi2r2(float f_pi, float fs)
{
    float d = -1 * cos((2 * PI * f_pi) / fs);
    float r2 = acos(-d);
    return r2;
}
float f_b2r1(float f_b, float fs)
{
    float tmp = tan(PI * f_b / fs);
    float c = (tmp - 1) / (tmp + 1);
    float r1 = acos(-c);
    return r1;
}

static void tvap_set_f_pi(t_tvap_tilde *x, t_floatarg f_pi)
{
    x->x_f_pi = f_pi;
}

static void tvap_set_f_b(t_tvap_tilde *x, t_floatarg f_b)
{
    x->x_f_b = f_b;
}

static t_int *tvap_tilde_perform(t_int *w)
{   
    float _z1, _z2, _r1, _r2, _cr1, _sr1, _cr2, _sr2;
    float tmp[3];

    t_sample *in   = (t_sample *)(w[1]);
    t_sample *f_pi = (t_sample *)(w[2]);
    t_sample *f_b  = (t_sample *)(w[3]);
    t_sample *out  = (t_sample *)(w[4]);
    t_tvap_ctl *ctl  = (t_tvap_ctl *)(w[5]);
    int n = (int)(w[6]);

    _z1 = ctl->z1;
    _z2 = ctl->z2;
    while (n--)
    {
    	float   _x_n = *(in++), 
                _f_pi_n = *(f_pi++), 
                _f_b_n = *(f_b++),
                _y_n;
        _r1 = f_b2r1(_f_b_n, ctl->fs);
        _r2 = f_pi2r2(_f_pi_n, ctl->fs);

        _cr1 = cos(_r1);
        _sr1 = sin(_r1);
        _cr2 = cos(_r2);
        _sr2 = sin(_r2);

        //tmp[0] = _x_n;
        tmp[1] = _cr2 * _z1 - _sr2 * _z2;
        //tmp[2] = _sr2 * _z1 + _cr2 * _z2;

        tmp[0] = _cr1 * _x_n - _sr1 * tmp[1];
        tmp[1] = _sr1 * _x_n + _cr1 * tmp[1];
        tmp[2] = _sr2 * _z1 + _cr2 * _z2;

        _z1 = tmp[1];
        _z2 = tmp[2];

        *out++ = tmp[0];
    }

    ctl->z1 = _z1;
    ctl->z2 = _z2;
    return (w+7);
}

static void tvap_tilde_dsp(t_tvap_tilde *x, t_signal **sp)
{
    x->x_sr = sp[0]->s_sr;
    x->x_ctl->fs = sp[0]->s_sr; // could be les redundant somehow...
    dsp_add(tvap_tilde_perform, 6, 
        sp[0]->s_vec, // input
        sp[1]->s_vec, // f_pi
        sp[2]->s_vec, // f_b
        sp[3]->s_vec, // output
        x->x_ctl,     // a struct to store the z1 and z2 history
        sp[0]->s_n);
}

void tvap_tilde_setup(void)
{
    tvap_tilde_class = class_new(gensym("tvap~"), (t_newmethod)tvap_tilde_new, 0,
    	sizeof(t_tvap_tilde), 0, A_DEFFLOAT, 0);
	    /* this is magic to declare that the leftmost, "main" inlet
	    takes signals; other signal inlets are done differently... */
    CLASS_MAINSIGNALIN(tvap_tilde_class, t_tvap_tilde, x_f);
    	/* here we tell Pd about the "dsp" method, which is called back
	when DSP is turned on. */
    class_addmethod(tvap_tilde_class, (t_method)tvap_tilde_dsp, 
        gensym("dsp"), A_CANT, 0);
}
