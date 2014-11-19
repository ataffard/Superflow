#pragma once

#include <math.h>

// Root Packages
#include "TLorentzVector.h"

namespace PhysicsTools {
    //Efficiency using binomial error
    void binomialError(double Num, double Den, double& Eff, double& EffErr);

    //
    //Physics quantities
    //

    //Transverse mass
    double mT(TLorentzVector* _l, TLorentzVector _nu);

    //Transverse mass WW
    double mTWW(TLorentzVector _ll, TLorentzVector _nu, bool MvvTrue = true);

    // Signed d0 wrt to jet axis 
    double signedD0(double d0, double sigmaD0, TLorentzVector _p, TLorentzVector _j);

    // pT relative to jet axis 
    double ptRel(TLorentzVector j, TLorentzVector p);


    /*
      References: mCT.0802.2879v3 & mCT.0910.1584v2
      */
    //mCT
    double mCT(TLorentzVector v1, TLorentzVector v2);

    //mCTperp
    double mCTperp(TLorentzVector lep0, TLorentzVector lep1, TLorentzVector met);

    //mCTpara
    double mCTpara(TLorentzVector lep0, TLorentzVector lep1, TLorentzVector met);

    //mll collinear approx
    double mZTauTau(TLorentzVector l0, TLorentzVector l1, TLorentzVector met);

    //Acoplananrity
    double acoplanarity(TLorentzVector v0, TLorentzVector v1);

    //mColl - LFV
    double mColl(TLorentzVector* lep0, TLorentzVector* lep1, TLorentzVector met);

    //Sphericity
    float Sphericity(std::vector<float> &pxVector,
        std::vector<float> &pyVector,
        std::vector<float> &pzVector,
        bool IsTransverseSphericity);
}