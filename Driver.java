
import FoilSim.Foil;
import Acme.MainFrame;

public class Driver {
    public static void main(String[] args){
        Foil mFoil = new Foil();
        MainFrame mMainFrame = new MainFrame(mFoil, 750, 500);
    }
}