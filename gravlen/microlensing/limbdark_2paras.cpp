double VBBinaryLensing::BinaryMagDark(double a, double q, double y1, double y2, double RSv, double a1, double Tol) {
	double Mag = -1.0, Magold = 0., Tolv = Tol;
	double tc, lb, rb, lc, rc, cb, cc, r2, cr2, scr2;
	int c = 0, flag;
	double currerr, maxerr;
	annulus *first, *scan, *scan2;
	int nannold, totNPS = 1;
	_sols *Images;

	y_1 = y1;
	y_2 = y2;
	while ((Mag < 0.9) && (c < 3)) {

		first = new annulus;
		first->bin = 0.;
		first->cum = 0.;// cumulative function F(r)
		if (Mag0 > 0.5) {
			first->Mag = Mag0;
			first->nim = nim0;
		}
		else {
			first->Mag = BinaryMag0(a, q, y_1, y_2, &Images);
			first->nim = Images->length;
			delete Images;
		}
		first->f = 3 / (3 - a1);// r = 0, so using point source mag computation above, Images is a _sols object, its length = nim -- number of image tracks
		first->err = 0;
		first->prev = 0;


		first->next = new annulus;
		scan = first->next;
		scan->prev = first;
		scan->next = 0;
		scan->bin = 1.;
		scan->cum = 1.;
		scan->Mag = BinaryMag(a, q, y_1, y_2, RSv, Tolv, &Images);
		totNPS += NPS;
		scan->nim = Images->length;
		delete Images;
		scan->f = first->f * (1 - a1); // r = 1, so using FS mag computation above
		if (scan->nim == scan->prev->nim) {
			scan->err = fabs((scan->Mag - scan->prev->Mag) * (scan->prev->f - scan->f) / 4);
		}
		else {
			scan->err = fabs((scan->Mag) * (scan->prev->f - scan->f) / 4);
		}

		Magold = Mag = scan->Mag;
		//			scan->err+=scan->Mag*Tolv*0.25; //Impose calculation of intermediate annulus at mag>4. Why?
		currerr = scan->err;
		flag = 0;
		nannuli = nannold = 1;
		while (((flag < nannold + 5) && (currerr > Tolv) && (currerr > RelTol * Mag)) || (nannuli < minannuli)) {
			// update/find the annuli with the maximum error
			maxerr = 0;
			for (scan2 = first->next; scan2; scan2 = scan2->next) {
#ifdef _PRINT_ERRORS_DARK
				printf("\n%d %lf %le | %lf %le", nannuli, scan2->Mag, scan2->err, Mag, currerr);
#endif
				if (scan2->err > maxerr) {
					maxerr = scan2->err;
					scan = scan2;
				}
			}
			nannuli++;
			Magold = Mag;
			Mag -= (scan->Mag * scan->bin * scan->bin - scan->prev->Mag * scan->prev->bin * scan->prev->bin) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
			currerr -= scan->err;
			rb = scan->bin; // current bin rb, right bin
			rc = scan->cum; // right cumulative function
			lb = scan->prev->bin; // prev bin lb, left bin
			lc = scan->prev->cum; // left cumulative function
			tc = (lc + rc) / 2; // Mid Cumulative Function as target
			do {
				cb = rb + (tc - rc) * (rb - lb) / (rc - lc); // center bin, bin actually is the corresponding radius in range (0,1)
				r2 = cb * cb;
				cr2 = 1 - r2;
				scr2 = sqrt(cr2);
				cc = (3 * r2 * (1 - a1) - 2 * a1 * (scr2 * cr2 - 1)) / (3 - a1); // cumulative function F(r) from equ(45) at new r = cb, integral of f(r)
				if (cc > tc) {
					rb = cb;
					rc = cc;
				}
				else {
					lb = cb;
					lc = cc;
				}
			} while (fabs(cc - tc) > 1.e-5); // loop until equal partition the cumulative function
			scan->prev->next = new annulus;
			scan->prev->next->prev = scan->prev;
			scan->prev = scan->prev->next;
			scan->prev->next = scan;
			scan->prev->bin = cb;
			scan->prev->cum = cc;
			scan->prev->f = first->f * (1 - a1 * (1 - scr2));//Bozza 2010 42
			scan->prev->Mag = BinaryMag(a, q, y_1, y_2, RSv * cb, Tolv, &Images);
			totNPS += NPS;
			scan->prev->nim = Images->length;
			if (scan->prev->prev->nim == scan->prev->nim) {
				scan->prev->err = fabs((scan->prev->Mag - scan->prev->prev->Mag) * (scan->prev->prev->f - scan->prev->f) * (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin) / 4);
			}
			else {
				scan->prev->err = fabs((scan->prev->bin * scan->prev->bin * scan->prev->Mag - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->Mag) * (scan->prev->prev->f - scan->prev->f) / 4);
			}
			if (scan->nim == scan->prev->nim) {
				scan->err = fabs((scan->Mag - scan->prev->Mag) * (scan->prev->f - scan->f) * (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin) / 4);
			}
			else {
				scan->err = fabs((scan->bin * scan->bin * scan->Mag - scan->prev->bin * scan->prev->bin * scan->prev->Mag) * (scan->prev->f - scan->f) / 4);
			}
			rb = (scan->Mag + scan->prev->prev->Mag - 2 * scan->prev->Mag);
			scan->prev->err += fabs(rb * (scan->prev->prev->f - scan->prev->f) * (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin));
			scan->err += fabs(rb * (scan->prev->f - scan->f) * (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin));
#ifdef _PRINT_ERRORS_DARK
			printf("\n%d", Images->length);
#endif
			delete Images;

			Mag += (scan->bin * scan->bin * scan->Mag - cb * cb * scan->prev->Mag) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
			Mag += (cb * cb * scan->prev->Mag - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->Mag) * (scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin);
			currerr += scan->err + scan->prev->err;

			if (fabs(Magold - Mag) * 2 < Tolv) {
				flag++;
			}
			else {
				flag = 0;
				nannold = nannuli;
			}

		}

		if (multidark) {
			annlist = first;
		} else {
			while (first) {
				scan = first->next;
				delete first;
				first = scan;
			}
		}

		Tolv /= 10;
		c++;
	}
	NPS = totNPS;
	therr = currerr;

	return Mag;
}
