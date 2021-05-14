encode(m)=fromdigits(Vec(Vecsmall(m)),128);
[n,c] = readvec("input.txt");


decode(m) = Strchr(digits(m, 128));

pollard_p_1(nb, borne) = {

  my(a);

  a = Mod(2 + random(nb-2), nb);    \\ nombre aléatoire entre 2 et nb, 0 et 1 inutiles.

  pollard_p(nb, borne, a);
}


pollard_p(nb, borne, a) = {

  my(aM, m, pgcd);

  aM = Mod(a, nb);

  forprime(p = 2, borne,
    m = p ^ logint(nb, p);
    aM = aM ^ m ;
    pgcd = gcd(lift(aM) -1, nb);

    if(pgcd == nb,                                  \\ Si pgcd = nb, alors "borne" trop grande.
      return( pollard_p(nb, p - 1, a ^ m) ),        \\ a ^ m car p − 1 = kl où k est "borne"-friable et "borne" < l < "borne" semi-friable
                                                    \\ donc on calcul a ^ m pour tout "borne" < m < " borne' ".
    \\ else
      if(pgcd != 1,                                 \\ Si pgcd = 1, alors "borne" trop petite.
        return(pgcd)
      )
    );
  );
}


Rabin(c, n, p, q) = {

  my(chiffre, jacobi, parite, m);

  chiffre = c[1];
  jacobi = c[2];
  parite = c[3];

  mp = Mod(chiffre, p)^( (p+1)/4 );
  mq = Mod(chiffre, q)^( (q+1)/4 );

	if(kronecker(lift(mp), p) == -1,
		mp = -mp
	);

	if(kronecker(lift(mq), q) == -1,
		mq = -mq
	);

	if(jacobi == 1,
		m = chinese(mp, mq);
		if(lift(m) % 2 == parite,
		  m,
    \\ else
			chinese(-mp, -mq)
		),
  \\ else
		m = chinese(-mp, mq);

		if(lift(m) % 2 == parite,
			m,
    \\ else
			chinese(mp, -mq)
		)
	);
};



main(nb, chiffre) = {

  my(p, q);

  p = pollard_p_1(nb, 10000);
  q = nb/p;

  decode( lift( Rabin(chiffre, nb, p, q) ) );
}

print(main(n, c));
