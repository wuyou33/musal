#ifndef MusAL_ExplicitButcherTabs_h
#define MusAL_ExplicitButcherTabs_h


/***
* Provides Butcher tableaux for explicit RK methods
*/
class ExplicitButcherTabs{

public:
	/* Default constructor */
    ExplicitButcherTabs();

    /* Destructor */
	~ExplicitButcherTabs();


	/**
	* Get methods
	*/
	/* Get ERK4 data */
	inline int gets_erk_4(){
		return s_erk_4;
	}

	inline double* getA_erk_4(){
		return A_erk_4;
	}

	inline double* getb_erk_4(){
		return b_erk_4;
	}

	inline double* getc_erk_4(){
		return c_erk_4;
	}

	/* Get explicit trapezoidal data */
	inline int gets_erk_trap(){
		return s_erk_trap;
	}

	inline double* getA_erk_trap(){
		return A_erk_trap;
	}

	inline double* getb_erk_trap(){
		return b_erk_trap;
	}

	inline double* getc_erk_trap(){
		return c_erk_trap;
	}

	/* Get Runge data */
	inline int gets_erk_rung(){
		return s_erk_rung;
	}

	inline double* getA_erk_rung(){
		return A_erk_rung;
	}

	inline double* getb_erk_rung(){
		return b_erk_rung;
	}

	inline double* getc_erk_rung(){
		return c_erk_rung;
	}


private:
    /* Explicit Trapezoidal */
    const static int s_erk_trap = 2;

    double A_erk_trap[4];
    double b_erk_trap[2];
    double c_erk_trap[2];

    /* Runge */
    const static int s_erk_rung = 2;

    double A_erk_rung[4];
    double b_erk_rung[2];
    double c_erk_rung[2];
    
	/* ERK4 */
	const static int s_erk_4 = 4;

	double A_erk_4[16];
	double b_erk_4[4];
	double c_erk_4[4];

};


#endif