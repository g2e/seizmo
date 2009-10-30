function [bad,bad2]=badsacpz(sacpz)

networks=fieldnames(sacpz);
idx=0; idx2=0;
for m=1:numel(networks)
    disp(networks{m})
    for n=1:numel(sacpz.(networks{m}))
        try
            cplxpair(sacpz.(networks{m})(n).p);
            cplxpair(sacpz.(networks{m})(n).z);
        catch
            idx=idx+1;
            bad(idx).net=networks{m};
            bad(idx).idx=n;
            bad(idx).sacpz=sacpz.(networks{m})(n);
            %{
            % need to decide which is the typo
            p=[-0.1053+7.864i
               -10.533-7.864i];
            % need to figure out what to do with -628.3 +     0.0016i
            p=[ -8.796 +      8.974i
                -8.796 -      8.974i
                -139.8 +      612.6i
                -391.7 +      491.2i
                -566.1 +      272.6i
                -628.3 +     0.0016i
                -566.1 -      272.6i
                -391.8 -      491.2i
                -139.8 -      612.6i];
            % just conjugate
            p=[-0.9211 +       0.94i
               -0.9211 +       0.94i];
            % just conjugate
            p=[   -981 +       1009i
                  -981 -       1009i
                 -3290 -       1263i
                 -3290 -       1263i];
            % looks like drop a line and add conj
            p=[-0.0123 +     0.0123i
               -0.0123 -     0.0123i
               -0.0123 -     0.0123i
               -19.588 +     24.562i];
            %
            p=[];
            %
            p=[];
            %
            p=[];
            %
            p=[];
            %
            p=[];
            %
            p=[];
            %
            p=[];
            %
            p=[];
            %}
        end
        try
            rpi=sacpz.(networks{m})(n).p==real(sacpz.(networks{m})(n).p);
            zpi=sacpz.(networks{m})(n).p==0;
            rp=sort(sacpz.(networks{m})(n).p(rpi & ~zpi));
            if(~isequal(rp,unique(rp)))
                error('badsacpz:doublepole','repeat non-zero poles!');
            end
            rzi=sacpz.(networks{m})(n).z==real(sacpz.(networks{m})(n).z);
            zzi=sacpz.(networks{m})(n).z==0;
            rz=sort(sacpz.(networks{m})(n).z(rzi & ~zzi));
            if(~isequal(rz,unique(rz)))
                error('badsacpz:doublezero','repeat non-zero zeros!');
            end
        catch
            idx2=idx2+1;
            bad2(idx2).net=networks{m};
            bad2(idx2).idx=n;
            bad2(idx2).sacpz=sacpz.(networks{m})(n);
        end
    end
end