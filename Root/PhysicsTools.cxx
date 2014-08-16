#include "Superflow/PhysicsTools.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include "TGraph.h"
#include "TMatrixD.h"
#include "TVectorD.h"

namespace PhysicsTools {
    //-----------------------------------------------------------------------------
    void binomialError(double Num, double Den, double& Eff, double& EffErr)
    {
        //Compute Eff=N1/N2 where N1=pass, N2=total (pass+fail)
        Eff = 0;
        EffErr = 0;
        Eff = (Den > 0. ? Num / Den : 0.);
        EffErr = (Den > 0. ? sqrt((Eff - pow(Eff, 2)) / Den) : 0.);
    }


    //-----------------------------------------------------------------------------
    // Transverse mass 
    double mT(TLorentzVector* _l, TLorentzVector _nu)
    {
        return sqrt(2 * _l->Pt()*_nu.Pt() * (1 - cos(_l->DeltaPhi(_nu))));
    }

    //-----------------------------------------------------------------------------
    // Transverse mass WW
    double mTWW(TLorentzVector _ll, TLorentzVector _nu, bool MvvTrue)
    {
        double dphi = acos(cos(_ll.Phi() - _nu.Phi()));
        double mll = _ll.M();
        double mvv = 0;
        if (!MvvTrue) mvv = mll;

        double mT = 0;
        mT = sqrt(pow(mll, 2) + pow(mvv, 2)
            + 2 * (sqrt(pow(_ll.Pt(), 2) + pow(mll, 2)) * sqrt(pow(_nu.Pt(), 2) + pow(mvv, 2))
            - _ll.Pt()*_nu.Pt()*cos(dphi)));

        return mT;
    }
    //-----------------------------------------------------------------------------
    double mCT(TLorentzVector v1, TLorentzVector v2)
    {
        double mct = (v1.Et() + v2.Et())*(v1.Et() + v2.Et()) - (v1 - v2).Perp2();
        return (mct < 0) ? 0.0 : sqrt(fabs(mct)); // AT
    };

    //-----------------------------------------------------------------------------
    double mCTperp(TLorentzVector lep0, TLorentzVector lep1, TLorentzVector met)
    {
        TVector3 v13D = lep0.Vect();
        TVector3 v23D = lep1.Vect();
        TVector3 u3D = met.Vect();
        u3D.SetZ(0);
        u3D = u3D.Unit();

        // Calculate p[1,2]Perp
        TVector3 p1Perp = u3D.Cross(v13D.Cross(u3D));
        TVector3 p2Perp = u3D.Cross(v23D.Cross(u3D));

        double ET1Perp = sqrt(p1Perp.Perp2() + lep0.M2());
        double ET2Perp = sqrt(p2Perp.Perp2() + lep1.M2());

        double mCTSQ = pow(ET1Perp + ET2Perp, 2) - (p1Perp - p2Perp).Perp2();

        return ((mCTSQ > 0.) ? sqrt(mCTSQ) : 0.);


    };

    //-----------------------------------------------------------------------------
    double mCTpara(TLorentzVector lep0, TLorentzVector lep1, TLorentzVector met)
    {
        TVector3 v13D = lep0.Vect();
        TVector3 v23D = lep1.Vect();
        TVector3 u3D = met.Vect();
        u3D.SetZ(0);
        u3D = u3D.Unit();

        // Calculate p[1,2]Para
        TVector3 p1Para = v13D.Dot(u3D)*u3D;
        TVector3 p2Para = v23D.Dot(u3D)*u3D;

        double ET1Para = sqrt(p1Para.Perp2() + lep0.M2());
        double ET2Para = sqrt(p2Para.Perp2() + lep1.M2());

        double mCTSQ = pow(ET1Para + ET2Para, 2) - (p1Para - p2Para).Perp2();

        return ((mCTSQ > 0.) ? sqrt(mCTSQ) : 0.);

    };

    //-----------------------------------------------------------------------------
    double mColl(TLorentzVector* lep0, TLorentzVector* lep1, TLorentzVector met)
    {
        double dEta = lep0->Eta() - lep1->Eta();
        double dPhi = lep0->DeltaPhi(*lep1);

        return sqrt(2 * lep0->Pt()*(lep1->Pt() + met.Et())*(cosh(dEta) - cos(dPhi)));
    }



    //-----------------------------------------------------------------------------
    // d0 signed wrt to jet direction
    double signedD0(double d0, double sigmaD0,
        TLorentzVector _p, TLorentzVector _j)
    {
        double sD0 = -999;
        double qd0 = d0 / fabs(d0);
        double m_sPhi = _p.Phi() + qd0 * TMath::Pi() / 2.;
        double dPhi = m_sPhi - _j.Phi();
        double signIP = fabs(cos(dPhi)) / (cos(dPhi) + 1e-32) * fabs(d0);
        sD0 = signIP / sigmaD0;

        return sD0;

    }
    /*--------------------------------------------------------------------------------*/
    // pTshik relshii to jet axis
    double ptRel(TLorentzVector j, TLorentzVector p)
    {
        TVector3 jet(j.Px(), j.Py(), j.Pz());
        TVector3 part(p.Px(), p.Py(), p.Pz());
        return part.Perp(jet);
    }

    /*--------------------------------------------------------------------------------*/
    // mtt collinear approximation
    //https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Summer
    //https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG4/software/common/CommonAnalysisUtils/trunk/Root/Algorithms.cxx
    double mZTauTau(TLorentzVector l0, TLorentzVector l1, TLorentzVector met)
    {
        double px0(l0.Px()), py0(l0.Py());
        double px1(l1.Px()), py1(l1.Py());
        double pxm(met.Px()), pym(met.Py());
        double num(px0*py1 - py0*px1);
        double den1(py1*pxm - px1*pym + px0*py1 - py0*px1);
        double den2(px0*pym - py0*pxm + px0*py1 - py0*px1);
        double x1 = (den1 != 0.0 ? (num / den1) : 0.0);
        double x2 = (den2 != 0.0 ? (num / den2) : 0.0);
        // not guaranteed that this configuration is kinematically possible
        return (x1*x2 > 0.0 ? (l0 + l1).M() / sqrt(x1*x2) : -1.0);
    }

    /*--------------------------------------------------------------------------------*/
    double acoplanarity(TLorentzVector v0, TLorentzVector v1)
    {
        return v1.Phi() - v0.Phi() - TMath::Pi();
    }


    //-----------------------------------------------------------------------------
    //Take from https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/SUSYD3PDSnippets#Transverse_Sphericity_ST
    // Fill the 3 vectors with px,py,pz of the different particles entering sphericity computation.
    float Sphericity(std::vector<float> &pxVector,
        std::vector<float> &pyVector,
        std::vector<float> &pzVector,
        bool IsTransverseSphericity)
    {
        //!matrix elements are set to 0
        // Defines the Sphericity Matrix
        TMatrixD sphericityMatrix(3, 3);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                sphericityMatrix[i][j] = 0.;
        Double_t d_sph_p = 0.;

        //////////////////////////
        //! Calculates Sphericity Matrix S_{ij} = \sum_k p_k^j p_k^j  / \sum_k |p_k|^2
        Double_t d_sph_Pxyz[3];
        for (int i = 0; i < 3; i++) d_sph_Pxyz[i] = 0.;
        for (unsigned int indexI = 0; indexI < pxVector.size(); indexI++) {
            for (int i = 0; i < 3; i++) d_sph_Pxyz[i] = 0.;
            d_sph_Pxyz[0] = pxVector[indexI];
            d_sph_Pxyz[1] = pyVector[indexI];

            //! only if the 3D Sph is calculated, one needs pz != 0
            if (!IsTransverseSphericity) d_sph_Pxyz[2] = pzVector[indexI];

            //! The Sphericity Matrix is calculated. If pz = 0 it
            //! is the needed 2x2 Matrix
            for (int iMatrixI = 0; iMatrixI < 3; iMatrixI++)
                for (int iMatrixJ = 0; iMatrixJ < 3; iMatrixJ++)
                    sphericityMatrix[iMatrixI][iMatrixJ] += d_sph_Pxyz[iMatrixI] * d_sph_Pxyz[iMatrixJ];

            //! calculates \sum_k |p_k|^2
            d_sph_p += d_sph_Pxyz[0] * d_sph_Pxyz[0] +
                d_sph_Pxyz[1] * d_sph_Pxyz[1] +
                d_sph_Pxyz[2] * d_sph_Pxyz[2];
        }

        //!  Normalizes S_{ij}
        if (d_sph_p != 0.) {
            for (int iMatrixI = 0; iMatrixI < 3; iMatrixI++) {
                for (int iMatrixJ = 0; iMatrixJ < 3; iMatrixJ++) {
                    sphericityMatrix[iMatrixI][iMatrixJ] =
                        sphericityMatrix[iMatrixI][iMatrixJ] / d_sph_p;
                }
            }
            //! if there are no values available, it crashes.
        }
        else {
            // Cleaning
            return -99.;
        }

        //! Calculate the EigenValues
        TVectorD eigenValues;
        const TMatrixD eigenVectoren = sphericityMatrix.EigenVectors(eigenValues);

        //! The EigenValues have to be sorted: Lambda1 > Lambda2 > Lambda3
        Int_t eigenLambda1 = 0, eigenLambda2 = 1, eigenLambda3 = 2;
        // from the babar sphericity code...
        double Emax = eigenValues[0];
        double Emin = eigenValues[0];
        for (int i = 0; i <= 2; ++i) {
            if (eigenValues[i] > Emax) {
                Emax = eigenValues[i];
                eigenLambda1 = i;
            }
            if (eigenValues[i] < Emin) {
                Emin = eigenValues[i];
                eigenLambda3 = i;
            }
        }
        eigenLambda2 = 3 - eigenLambda3 - eigenLambda1;

        //! Calculates the Sphericity with
        //! S_T = 2 \lambda_2/(\lambda_1 + \lambda_2)
        //! TDR Vol II, May 1999, Page 820
        double sphericity = -99.;
        if (IsTransverseSphericity)
            sphericity = 2.*eigenValues[eigenLambda2] /
            (eigenValues[eigenLambda1] + eigenValues[eigenLambda2]);
        else
            sphericity = (3. / 2.)*
            (eigenValues[eigenLambda2] + eigenValues[eigenLambda3]);

        return sphericity;
    }
}