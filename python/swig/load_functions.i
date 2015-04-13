%inline %{
  namespace akantu {

    static void lin_load(double * position, double * load, 
		       __attribute__ ((unused)) Real * normal,
		       __attribute__ ((unused)) UInt surface_id) {
      
      memset(load,0,sizeof(Real)*3);
      if (position[0]<=10){
	load[1]= -6000;
      }
    }
  }
  %}
