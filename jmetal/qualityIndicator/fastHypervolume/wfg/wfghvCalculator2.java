package jmetal.qualityIndicator.fastHypervolume.wfg;

import java.io.IOException;

import jmetal.core.Problem;
import jmetal.core.*;

public class wfghvCalculator2 {
	public jmetal.qualityIndicator.util.MetricsUtil utils_;
	SolutionSet pf_ = null;
	double [][] pfMatrix_ = null;
	int fun;
	public wfghvCalculator2(SolutionSet paretoFront) {
	    pf_ = paretoFront;
	    pfMatrix_ = null;  
	    utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	  } // Constructor 
	
	public wfghvCalculator2(SolutionSet paretoFront,int fun) {
	    pf_ = paretoFront;
	    pfMatrix_ = null;
	    this.fun = fun;
	    utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	  } // Constructor
	public static double hv2point(Solution Point1, Solution ref){
		  double x=ref.getObjective(0)-Point1.getObjective(0);
			for (int j=1;j<Point1.numberOfObjectives();j++){
				 x = x*(ref.getObjective(j)-Point1.getObjective(j));
			}
			return x;
		}
	public double calculatewfghv() throws IOException{

	    double hv;
	    int number = pf_.get(0).numberOfObjectives();
	    
	    Solution referencePoint1 = new Solution( number);
		
		for (int j=0;j<number;j++){
			if(fun == 6){
				referencePoint1.setObjective(j,0.5);
			}else if(fun>6&fun<=9){
				referencePoint1.setObjective(j,1.0);
			}else if(fun>9&fun<=11){
				if(j!=number-1){
					referencePoint1.setObjective(number-1-j,Math.pow(Math.sqrt(2)/2, j));
				}else{
					referencePoint1.setObjective(0,referencePoint1.getObjective(1));
				}
			}else if(fun>12&fun<=21){
				referencePoint1.setObjective(j,2.0*(j+1));
			}
		}
		
		//NORMALIZATION
		for (int j=0;j<pf_.size();j++)
			for(int k=0;k<number;k++)
				pf_.get(j).setObjective(k, pf_.get(j).getObjective(k)/(1.1*referencePoint1.getObjective(k)) );
		SolutionSet invertedFront;
		//invertedFront = utils_.invertedFront(pf_,number);
		for (int j=0;j<pf_.size();j++)
			for(int k=0;k<number;k++)
				if(pf_.get(j).getObjective(k)>1.0){
					//pf_.remove(j);
					//j--;
					for(int s=0;s<number;s++){
						pf_.get(j).setObjective(s, 1.0);
					}
					break;
				}
		for (int j=0;j<number;j++)
			//referencePoint1.setObjective(j,2.0*j+10);
		    referencePoint1.setObjective(j,1.0);
		if (pf_.size() == 0){
			hv = 0.0;
		}else if(pf_.size() == 1){
			hv = hv2point(pf_.get(0), referencePoint1);
		}else if (pf_.size() == 2){
			double [] mid = new double [number];
			for (int j=0;j<number;j++){
				mid[j] = Math.max(pf_.get(0).getObjective(j), pf_.get(1).getObjective(j));
			}
			Solution midp = new Solution(number);
			for(int i=0;i<number;i++){
				midp.setObjective(i, mid[i]);
			}
			hv = hv2point(pf_.get(0),referencePoint1)+hv2point(pf_.get(1),referencePoint1)-hv2point(midp,referencePoint1);
			
		}else if (pf_.size()==3){
			double [] w01 = new double [number];
			double [] w02 = new double [number];
			double [] w12 = new double [number];
			double [] w012 = new double [number];
			for (int j=0;j<number;j++){
				w01[j] = Math.max(pf_.get(0).getObjective(j), pf_.get(1).getObjective(j));
				w02[j] = Math.max(pf_.get(0).getObjective(j), pf_.get(2).getObjective(j));
				w12[j] = Math.max(pf_.get(1).getObjective(j), pf_.get(2).getObjective(j));
			}
			for (int j=0;j<number;j++){
				w012[j] = Math.max(w02[j], pf_.get(1).getObjective(j));
			}
			Solution p01 = new Solution(number);Solution p02 = new Solution(number);Solution p12 = new Solution(number);Solution p012 = new Solution(number);
			for(int i=0;i<number;i++){
				p01.setObjective(i, w01[i]);
				p02.setObjective(i, w02[i]);
				p12.setObjective(i, w12[i]);
				p012.setObjective(i, w012[i]);
			}
			hv = hv2point(pf_.get(0),referencePoint1)+hv2point(pf_.get(1),referencePoint1)+hv2point(pf_.get(2),referencePoint1)
					-hv2point(p01,referencePoint1)-hv2point(p02,referencePoint1)-hv2point(p12,referencePoint1)+hv2point(p012,referencePoint1);
		}else{
			WFGHV1 wfghv = new WFGHV1(number, pf_.size()) ;
			Front1 front = new Front1(pf_.size(), number, pf_);
		    hv = wfghv.getHV(front,referencePoint1);
		}
		return hv;
  } // CalculateHypervolume
}
