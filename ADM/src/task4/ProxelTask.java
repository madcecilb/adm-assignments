package task4;


public class ProxelTask {
	
	public class Proxel {
	    int     id;                  /* unique proxel id for searching    */
	    int     s;                   /* discrete state of SPN             */
	    int     tau1k;               /* first supplementary variable      */
	    int     tau2k;               /* second supplementary variable     */
	    double  val;                 /* proxel probability                */
	    Proxel first, second, third, healthy;         /* pointers to child proxels in tree */
	}
	
	public enum State {
	    Stage1, Stage2 , Stage3, Healthy
	}
	
	private static final double MINPROB = 1.0e-12;
	private static final double DELTA = 0.005;
	private static final int ENDTIME = 10;
	
	double[][]  y = new double[3][];             /* vectors for storing solution      */
	double  tmax;         /* maximum simulation time           */
	int     TAUMAX;
	int     totcnt;                 /* counts total proxels processed    */
	int     maxccp;                 /* counts max # concurrent proxels   */
	int     ccpcnt;                 /* counts concurrent proxels         */
	Proxel[] root = new Proxel[2];                /* trees for organising proxels      */
	Proxel[] firstfree;      /* linked list of free proxels       */
	double  eerror      = 0;         /* accumulated error                 */
	int     sw         = 0;         /* switch for old and new time steps */
	int len;
	double dt;
	
	
	/********************************************************/
	/*	distribution functions			                    */
	/*	instantaneous rate functions			            */
	/********************************************************/
	 

	/* returns weibull IRF */
	public double weibullhrf(Double x, Double alpha, Double beta, Double x0) {
	    return beta/alpha*Math.pow((x-x0)/alpha,beta-1);
	}


	/* returns deterministic IRF */
	double dethrf(double x, double d) {
	    double y;

	    if (Math.abs(x - d) < dt/2)
	        y = 1.0/dt;
	    else
	        y = 0.0;
	    return y;

	}


	/* returns uniform IRF */
	double unihrf(double x, double a, double b) {
	    double y;

	    if ((x >= a) && (x < b))
	        y = 1.0/(b-x);
	    else
	        y = 0.0;

	    return(y);
	}

	/* returns exponential IRF */
	double exphrf(double x, double l) {
	    return(l);
	}

	double normalpdf(double x, double m, double s){
	   double z = (x-m)/s;

	   return (Math.exp(-z * z/2)/(Math.sqrt(2*Math.PI)*s));
	}


	double logGamma(double x){
	   double coef[] = {76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.00120858003, -0.00000536382};
	   double step = 2.50662827465, fpf = 5.5, t, tmp, ser;
	   int i;

	   t = x - 1;
	   tmp = t + fpf;
	   tmp = (t + 0.5) * Math.log(tmp) -tmp;
	   ser = 1;
	   for(i = 1; i <= 6; i++){
	      t = t+1;
	      ser = ser + coef[i-1]/t;
	   }
	   return (tmp+Math.log(step * ser));
	}

	double gammaSeries(double x, double a){
	   int n, maxit = 100;
	   double eps = 0.0000003;
	   double sum = 1.0 / a, ap = a, gln = logGamma(a), del = sum;

	   for (n = 1; n <= maxit; n++){
	      ap++;
	      del = del * x / ap;
	      sum = sum + del;
	      if (Math.abs(del) < Math.abs(sum) * eps) break;
	   }
	   return (sum * Math.exp(-x + a * Math.log(x) - gln));
	}

	public double gammaCF(double x, double a){
	   int n, maxit = 100;
	   double eps = 0.0000003;
	   double gln = logGamma(a), g = 0, gold = 0, a0 = 1, a1 = x, b0 = 0, b1 = 1, fac = 1;
	   double an, ana, anf;

	   for(n = 1; n <= maxit; n++){
	      an = 1.0 * n;
	      ana = an - a;
	      a0 = (a1 + a0 * ana) * fac;
	      b0 = (b1 + b0 * ana) * fac;
	      anf = an * fac;
	      a1 = x * a0 + anf * a1;
	      b1 = x * b0 + anf * b1;
	      if (a1 != 0){
	         fac = 1.0 / a1;
	         g = b1 * fac;
	         if (Math.abs((g-gold)/g) < eps) 
	            break;
	         gold = g;
	      }
	   }
	   return (Math.exp(-x + a * Math.log(x) - gln) * g);
	}

	double gammacdf(double x, double a){
	   if (x <= 0) 
	      return 0;
	   else 
	   if (x < a + 1) 
	      return gammaSeries(x,a);
	   else 
	      return (1 - gammaCF(x, a));
	}

	double normalcdf(double x, double m, double s){
	   double z = (x - m) / s;

	   if (z >= 0)
	      return 0.5 + 0.5 * gammacdf(z * z / 2, 0.5);
	   else 
	      return (0.5 - 0.5 * gammacdf(z * z / 2, 0.5));
	}

	/* returns normal IRF */
	double normalhrf(double x, double m, double s){
	   return(normalpdf(x, m, s)/(1 - normalcdf(x, m, s)));
	}

	double lognormalpdf(double x, double a, double b){
	   double z = (Math.log(x) - a) / b;

	   return(Math.exp(- z * z / 2) / (x * Math.sqrt(2 * Math.PI) * b));
	}

	double lognormalcdf(double x, double a, double b){
	    double z = (Math.log(x) - a) / b;

	    if(x == 0)
	        return 0;
	    if (z >= 0) 
	        return(0.5 + 0.5 * gammacdf(z * z / 2, 0.5));
	    else 
	        return(0.5 - 0.5 * gammacdf(z * z / 2, 0.5));
	}

	/* returns lognormal IRF using mu & sigma */
	double lognormalhrf(double x, double a, double b){
	   if ((x == 0.0) || (x > 70000))
	      return(0.0);
	   else
	      return(lognormalpdf(x, a, b) / (1.0 - lognormalcdf(x, a, b)));
	}
	
	
	/********************************************************/
	/*	output functions			                        */
	/********************************************************/

	/* print all proxels in tree */
	void printtree(Proxel p) {
	    if (p == null)
	        return;
	    //String.format("%" + places + "." + decimals + "f", floatValue);
	    //printf("s %d t1 %d t2 %d val %lf \n",p.s,p.tau1k,p.tau2k,p.val);
	    System.out.println(p.s + "\t" + p.tau1k + "\t" + p.tau2k + "\t" + p.val);
	    printtree(p.first);
	    printtree(p.second);
	    printtree(p.third);
	    printtree(p.healthy);
	}

	/* print out complete solution */
	void plotsolution(int kmax) {
		System.out.println(" ");
		System.out.println(" ");
	    int k;

	    for(k=0; k<=kmax; k++)
	    	System.out.println(String.format("." + 5, k*dt, y[0][k], y[2][k]));
	        //printf("%7.5lf\t%7.5le\t%7.5le\n", k*dt, y[0][k], y[2][k]);

	}

	/* print out a proxel */
	void printproxel(Proxel c) {
		System.out.println(String.format("." + 5, c.s, c.tau1k,c.val));
	   // printf("processing %2d %2d %7.5le \n", c.s, c.tau1k, c.val);
	}
	
	/********************************************************/
	/*	proxel manipulation functions			            */
	/********************************************************/

	/* compute unique id from proxel state */
	int state2id(int s, int t1k, int t2k) {
	    return(TAUMAX*(TAUMAX*s+t1k)+t2k);
	}

	/* compute size of tree */
	int size(Proxel p) {
	    int s1, s2, s3, sh;

	    if (p == null)
	        return(0);
	    s1 = size(p.first);
	    s2 = size(p.second);
	    s3 = size(p.third);
	    sh = size(p.healthy);
	    return (s1+s2+s3+sh+1);
	}

	/* returns a proxel from the tree */
	public Proxel[] getproxel()
	{
	    Proxel temp;
	    Proxel old;
	    int LEFT = 0, RIGHT = 1;
	    int dir, cont = 1;
	    State dir1;
	    
	    if (root[1-sw] == null)
	        return(null);
	    temp = root[1-sw];
	    old  = temp;

	    /* move down the tree to a leaf */
	    while (cont == 1)
	    {
	        /* go first */
	        if ((temp.first != null) && (temp.second == null) && (temp.third == null) && (temp.healthy == null))
	        {
	            old  = temp;
	            temp = temp.first;
	            dir1  = State.Stage1;
	        }
	        /* go second */
	        else if ((temp.first == null) && (temp.second != null) && (temp.third == null) && (temp.healthy == null))
	        {
	            old  = temp;
	            temp = temp.second;
	            dir1  = State.Stage2;
	        }
	        /* go third */
	        else if ((temp.first == null) && (temp.second == null) && (temp.third != null) && (temp.healthy == null))
	        {
	            old  = temp;
	            temp = temp.third;
	            dir1  = State.Stage3;
	        }
	        /* go fourth */
	        else if ((temp.first == null) && (temp.second == null) && (temp.third == null) && (temp.healthy != null))
	        {
	            old  = temp;
	            temp = temp.healthy;
	            dir1  = State.Healthy;
	        }
	        /* choose right/left at random */
	        else if ((temp.right != NULL) && (temp->left != NULL))
	        {
	            if (rand() > RAND_MAX/2)
	            {
	                old  = temp;
	                temp = temp->left;
	                dir  = LEFT;
	            }
	            else
	            {
	                old  = temp;
	                temp = temp->right;
	                dir  = RIGHT;
	            }
	        }
	        else
	            cont = 0;
	    }
	    if (temp == root[1-sw])
	        root[1-sw] = NULL;
	    else
	    {
	        if (dir == RIGHT)
	            old->right = NULL;
	        else
	            old->left  = NULL;
	    }
	    old = firstfree;
	    firstfree = temp;
	    temp->right = old;
	    ccpcnt -= 1;
	    return(temp);
	}

	/* get a fresh proxel and copy data into it */
	proxel *insertproxel(int s, int tau1k, int tau2k, double val) {
	    proxel *temp;

	    /* create new proxel or grab one from free list */
	    if (firstfree == NULL)
	        temp = malloc(sizeof(proxel));
	    else {
	        temp = firstfree;
	        firstfree = firstfree->right;
	    }
	    /* copy values */
	    temp->id    = state2id(s, tau1k, tau2k);
	    temp->s     = s;
	    temp->tau1k = tau1k;
	    temp->tau2k = tau2k;
	    temp->val   = val;
	    ccpcnt     += 1;
	    if (maxccp < ccpcnt) {
	        maxccp = ccpcnt;
	        //printf("\n ccpcnt=%d",ccpcnt);
	    }
	    return(temp);
	}

	/* adds a new proxel to the tree */
	void addproxel(int s, int tau1k, int tau2k, double val) {
	    proxel *temp, *temp2;
	    int cont = 1,id;

	    /* Alarm! TAUMAX overstepped! */
	    if (tau1k >= TAUMAX) {
	        //  printf(">>> %3d %3d %3d %7.5le \n", s, tau1k, val, TAUMAX);
	        tau1k = TAUMAX - 1;
	    }


	    /* New tree, add root */
	    if (root[sw] == NULL) {
	        root[sw] = insertproxel(s,tau1k, tau2k, val);
	        root[sw]->left = NULL;
	        root[sw]->right = NULL;
	        return;
	    }

	    /* compute id of new proxel */
	    id = state2id(s,tau1k, tau2k);

	    /* Locate insertion point in tree */
	    temp = root[sw];
	    while (cont == 1) {
	        if ((temp->left != NULL) && (id < temp->id))
	            temp = temp->left;
	        else
	            if ((temp->right != NULL) && (id > temp->id))
	                temp = temp->right;
	            else
	                cont = 0;
	    }

	    /* Insert left leaf into tree */
	    if ((temp->left == NULL) && (id < temp->id)) {
	        temp2        = insertproxel(s, tau1k,tau2k, val);
	        temp->left   = temp2;
	        temp2->left  = NULL;
	        temp2->right = NULL;
	        return;
	    }

	    /* Insert right leaf into tree */
	    if ((temp->right == NULL) && (id > temp->id)) {
	        temp2        = insertproxel(s, tau1k,tau2k, val);
	        temp->right  = temp2;
	        temp2->left  = NULL;
	        temp2->right = NULL;
	        return;
	    }

	    /* Proxels have the same id, just add their vals */
	    if (id == temp->id) {
	        temp->val += val;
	        return;
	    }
	    printf("\n\n\n!!!!!! addproxel failed !!!!!\n\n\n");
	}
}



