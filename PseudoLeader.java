import comp34120.ex2.PlayerImpl;
import comp34120.ex2.PlayerType;
import comp34120.ex2.Record;
import java.rmi.NotBoundException;
import java.rmi.RemoteException;
import org.apache.commons.math3.linear.*;
import java.util.*;

/**
 * A pseudo leader. The members m_platformStub and m_type are declared
 * in the PlayerImpl, and feel free to use them. You may want to check
 * the implementation of the PlayerImpl. You will use m_platformStub to access
 * the platform by calling the remote method provided by it.
 * @author Xin
 */
final class PseudoLeader
	extends PlayerImpl
{
	Record[] records;
	double lambda = 0.952;
	double a,b;
	double[][] XT;
	double[] YT;
	double mape = 0;
	double rmse = 0;
	int d_count = 0;
	private int steps = 0;

	/**
	 * In the constructor, you need to call the constructor
	 * of PlayerImpl in the first line, so that you don't need to
	 * care about how to connect to the platform. You may want to throw
	 * the two exceptions declared in the prototype, or you may handle it
	 * by using "try {} catch {}". It's all up to you.
	 * @throws RemoteException
	 * @throws NotBoundException
	 */
	PseudoLeader()
		throws RemoteException, NotBoundException
	{
		/* The first parameter *MUST* be PlayerType.LEADER, you can change
		 * the second parameter, the name of the leader, such as "My Leader" */
		super(PlayerType.LEADER, "Pseudo Leader");
	}

	public static void main(final String[] p_args)	throws RemoteException, NotBoundException {
		new PseudoLeader();
	}

	/**
	 * You may want to delete this method if you don't want modify
	 * the original connection checking behavior of the platform.
	 * Actually I recommend you to delete this method from your own code
	 * @throws RemoteException If implemented, the RemoteException *MUST* be
	 * thrown by this method
	 */
	@Override
	public void checkConnection()
		throws RemoteException
	{
		super.checkConnection();
		//TO DO: delete the line above and put your own code here
	}

	/**
	 * You may want to delete this method if you don't want the platform
	 * to control the exit behavior of your leader class
	 * @throws RemoteException If implemented, the RemoteException *MUST* be
	 * thrown by this method
	 */
	@Override
	public void goodbye()
		throws RemoteException
	{
		super.goodbye();
		//TO DO: delete the line above and put your own exit code here
	}

	RealMatrix matrixP0;
	RealMatrix matrixTheta0;

	/**
	 * You may want to delete this method if you don't want to do any
	 * initialization
	 * @param p_steps Indicates how many steps will the simulation perform
	 * @throws RemoteException If implemented, the RemoteException *MUST* be
	 * thrown by this method
	 */
	@Override
	public void startSimulation(int p_steps) throws RemoteException	{
		try {
		records = new Record[130];
		steps = p_steps;
	/*	double[][] XT = {{1,3},
										{1,4},
										{1,5},
										{1,6},
										{1,7}};
		double[] YT = {2,3,3,4,6};*/
		XT = new double[130][2];
		YT = new double[130];
		for (int i = 1; i <= 99; i ++) {
				records[i] =	m_platformStub.query(m_type, i);
				XT[i-1][0] = 1;
				XT[i-1][1] = records[i].m_leaderPrice;
				YT[i-1] = records[i].m_followerPrice;
		}




		RealMatrix matrixXT = new Array2DRowRealMatrix(XT);
		RealMatrix matrixYT = new Array2DRowRealMatrix(YT);
		RealMatrix matrixXTT = matrixXT.transpose();

		RealMatrix theta = MatrixUtils.inverse(matrixXTT.multiply(matrixXT)).multiply(matrixXTT).multiply(matrixYT);
		double[][] thetaValues = theta.getData();
	/*	Uf = Rf(Ul)
		Uf = thetaValues[0]*

		Ul = (0.3b - 0.3a - 3)/(0.6b-2);

		Profit = (Ul - 1)*(2-Ul+0.3(thetaValues[0] + (thetaValues[1] * Ul)));
		Profit'= 0.6bUl - 2Ul + 3 + 0.3a -0.3b
		profit''=0.6b-2 */

		double[][] P0 = {
										{1/2000f, 0},
										{0, 1/2000f}
										};
		matrixP0 = new Array2DRowRealMatrix(P0);
			//m_platformStub.log(m_type, Arrays.deepToString(matrixP0.getData()));
		double[] theta0 = {1.5,0.140};
		matrixTheta0 = new Array2DRowRealMatrix(theta0);

		// go to 99 so proceedNewDay works for all days
		for (int i = 1; i <=99 ; i++) {
			RealMatrix temp = matrixP0;
			matrixP0 = calculatePtPlusOne(matrixP0, i);
			m_platformStub.log(m_type, Arrays.deepToString(matrixP0.getData()));
			RealMatrix LtPlusOne = calculateLtPlusOne(temp, i);
		//	m_platformStub.log(m_type, Arrays.deepToString(LtPlusOne.getData()));

			matrixTheta0 = updateParameterRLSA(matrixTheta0, calculateLtPlusOne(temp, i), i);
			//m_platformStub.log(m_type, Arrays.deepToString(matrixTheta0.getData()));
		}



	//	m_platformStub.log(m_type, Arrays.deepToString(theta.getData()));
//		m_platformStub.log(m_type, Arrays.deepToString(matrixTheta0.getData()));
	} catch (Exception e) {
		m_platformStub.log(m_type,e.toString());
	}
	}

	public RealMatrix calculateLtPlusOne(RealMatrix Pt, int day) throws RemoteException {
		double[] temp = {1, records[day].m_leaderPrice};
		RealMatrix phiXtPlusOne = new Array2DRowRealMatrix(temp);
		RealMatrix phiXtPlusOneTranspose = phiXtPlusOne.transpose();
	//	m_platformStub.log(m_type, Arrays.deepToString(phiXtPlusOne.getData()));
	//	m_platformStub.log(m_type, Arrays.deepToString(phiXtPlusOneTranspose.getData()));

		RealMatrix top = Pt.multiply(phiXtPlusOne);
		RealMatrix bottom = phiXtPlusOneTranspose.multiply(Pt).multiply(phiXtPlusOne).scalarAdd(lambda);
		//m_platformStub.log(m_type, Arrays.deepToString(top.getData()));
		//m_platformStub.log(m_type, Arrays.deepToString(bottom.getData()));

		return top.scalarMultiply(1f/(float)bottom.getTrace());
	}

	public RealMatrix calculatePtPlusOne(RealMatrix Pt, int day) {
		double[] temp = {1, records[day].m_leaderPrice};
		RealMatrix phiXtPlusOne = new Array2DRowRealMatrix(temp);
		RealMatrix phiXtPlusOneTranspose = phiXtPlusOne.transpose();

		RealMatrix top = Pt.multiply(phiXtPlusOne).multiply(phiXtPlusOneTranspose).multiply(Pt);
	//	m_platformStub.log(m_type, Arrays.deepToString(top.getData()));
		RealMatrix bottom = phiXtPlusOneTranspose.multiply(Pt).multiply(phiXtPlusOne).scalarAdd(lambda);
	//	m_platformStub.log(m_type, Arrays.deepToString(bottom.getData()));

		return Pt.subtract(top.scalarMultiply(1f/bottom.getTrace())).scalarMultiply(1f/(double)lambda);
	}

	public RealMatrix updateParameterRLSA(RealMatrix thetaT,RealMatrix LtPlusOne, int day) throws RemoteException{
		double[] temp = {1, records[day].m_leaderPrice};
		RealMatrix phiXtPlusOne = new Array2DRowRealMatrix(temp);
		RealMatrix phiXtPlusOneTranspose = phiXtPlusOne.transpose();
		double yTPlusOne = records[day].m_followerPrice;

		double predictionError = yTPlusOne - phiXtPlusOneTranspose.multiply(thetaT).getTrace();
			m_platformStub.log(m_type, "Prediction: " + Arrays.deepToString(phiXtPlusOneTranspose.multiply(thetaT).getData()));
	/*		m_platformStub.log(m_type, Double.toString(yTPlusOne));
			m_platformStub.log(m_type, Double.toString(predictionError)); */
		//m_platformStub.log(m_type, Arrays.deepToString(LtPlusOne.scalarMultiply(predictionError).getData()));


		return thetaT.add(LtPlusOne.scalarMultiply(predictionError));
	}

	/**
	 * You may want to delete this method if you don't want to do any
	 * finalization
	 * @throws RemoteException If implemented, the RemoteException *MUST* be
	 * thrown by this method
	 */
	@Override
	public void endSimulation()
		throws RemoteException
	{
		super.endSimulation();
		mape /= d_count;
		rmse /= d_count;
		rmse = Math.sqrt(rmse);
		m_platformStub.log(m_type, "MAPE: " + String.valueOf(mape));
		m_platformStub.log(m_type, "RMSE: " + String.valueOf(rmse));
		//TO DO: delete the line above and put your own finalization code here
		m_platformStub.log(m_type, "Total profit: " + calculateTotalProfit(steps));
	}

/*
	public void updateParameterRLSA(double Lt, double Pt, double[] XtPlusOne, double YtPlusOne, List<float> phiXt) {

		LtPlusOne = Pt*phiXt
		PtPlusOne
	} */

	public void updatePhiXt(List<Float> phiXt) {
		phiXt.add(records[phiXt.size()].m_leaderPrice);
	}

	/**
	 * To inform this instance to proceed to a new simulation day
	 * @param p_date The date of the new day
	 * @throws RemoteException This exception *MUST* be thrown by this method
	 */
	@Override
	public void proceedNewDay(int p_date)
		throws RemoteException
	{


				recursiveLeastSquared(p_date);
				//movingWindowApproach(p_date, 90);

				m_platformStub.log(m_type, "a: " + a);
				m_platformStub.log(m_type, "b: " + b);

				Double optimalPrice = ((0.3 * b)  - (0.3 * a) - 3)/((0.6 * b) - 2);

				m_platformStub.publishPrice(m_type, optimalPrice.floatValue());



		/*
		 * Check for new price
		 * Record l_newRecord = m_platformStub.query(m_type, p_date);
		 *
		 * Your own math model to compute the price here
		 * ...
		 * float l_newPrice = ....
		 *
		 * Submit your new price, and end your phase
		 * m_platformStub.publishPrice(m_type, l_newPrice);
		 */
	}

	public void movingWindowApproach(int p_date, int windowSize) throws RemoteException {
		p_date -= windowSize;
		p_date--;
		XT = new double[windowSize][2];
		YT = new double[windowSize];
		for (int i = 1; i <= windowSize; i++) {
				records[i] =	m_platformStub.query(m_type, i+p_date);
				XT[i-1][0] = 1;
				XT[i-1][1] = records[i].m_leaderPrice;
				YT[i-1] = records[i].m_followerPrice;
		}


		RealMatrix matrixXT = new Array2DRowRealMatrix(XT);
		RealMatrix matrixYT = new Array2DRowRealMatrix(YT);
		RealMatrix matrixXTT = matrixXT.transpose();

		RealMatrix theta = MatrixUtils.inverse(matrixXTT.multiply(matrixXT)).multiply(matrixXTT).multiply(matrixYT);
		double[][] thetaValues = theta.getData();
		a = thetaValues[0][0];
		b = thetaValues[1][0];

		double[] tempXt = {1, records[1].m_leaderPrice};
		RealMatrix phiXtPlusOne = new Array2DRowRealMatrix(tempXt).transpose();

		double pError = records[1].m_followerPrice - phiXtPlusOne.multiply(theta).getTrace();
		double pct_error = pError/records[1].m_followerPrice;
		mape += pct_error;
		rmse += (pError*pError);
		d_count++;
	}

	public void recursiveLeastSquared(int p_date) throws RemoteException {
		p_date--; //minus 1 to find to update from the previous days newly published data

		//m_platformStub.query(m_type, 1).m_followerPrice


		//get previous days data
		records[p_date] =	m_platformStub.query(m_type, p_date);
		//m_platformStub.log(m_type, "" + records[p_date].m_followerPrice);

		//update follower model with new data
		RealMatrix temp = matrixP0;
		matrixP0 = calculatePtPlusOne(matrixP0, p_date);
			m_platformStub.log(m_type, Arrays.deepToString(matrixP0.getData()));
		RealMatrix LtPlusOne = calculateLtPlusOne(temp, p_date);
		//	m_platformStub.log(m_type, Arrays.deepToString(LtPlusOne.getData()));

		matrixTheta0 = updateParameterRLSA(matrixTheta0, calculateLtPlusOne(temp, p_date), p_date);
		//m_platformStub.log(m_type, Arrays.deepToString(matrixTheta0.getData()));

		a = matrixTheta0.getEntry(0,0);
		b = matrixTheta0.getEntry(1,0);

		double[] tempXt = {1, records[1].m_leaderPrice};
		RealMatrix phiXtPlusOne = new Array2DRowRealMatrix(tempXt).transpose();

		double pError = records[1].m_followerPrice - phiXtPlusOne.multiply(matrixTheta0).getTrace();
		double pct_error = pError/records[1].m_followerPrice;
		mape += pct_error;
		rmse += (pError*pError);
		d_count++;
	}

	public double calculateTotalProfit(int p_steps) throws RemoteException
  {
    //for all dates, get records
    double total = 0.0;

    for(int i = 101; i <= (100 + p_steps); i++)
    {
      Record record = m_platformStub.query(m_type, i);
      total += (record.m_leaderPrice - 1) * ( 2 - record.m_leaderPrice + 0.3 * record.m_followerPrice);
    }

    return total;
  }
}
