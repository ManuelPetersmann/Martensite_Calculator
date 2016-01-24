function [ y1, y3, e1, e3  ] = ordered_eigenvalues( Ct )
% calculates the eigenvalues and eigenvectors and returns them in ascending
% order (of eigenvalues with correlating eigenvectors)

[W,D]=eig(Ct); % W...eigenvektoren e1, e3 von Ct korrespondierend zu eigenwerten
               % k1 und k3 von Ch, (D keine Anwendung hier, nur braucht
               % Befehl 2 Matrizen)

lambda_vorher = eig(Ct);
lambda_nachher = sort(lambda_vorher);% aufsteigende Sortierung,
% k1 <= k2 <= k3 muss gewährleistet sein!!!
y1 = lambda_nachher(1);
y2 = lambda_nachher(2);
y3 = lambda_nachher(3);

% if y1 > 1  || abs(1-y2) > 0.000001 
%     display('Eigenvalues do not satisfy conditions: y1>=0 and y2=0; it must be concluded that Ui and Uj are not Twin related')
%     not_existent = 1;
% else
%     not_existent = 0;
% end

%nachfolgend werden die korrespondierenden Eigenvektoren der Eigenwerte
%in der Matrix W wie die Eigenwerte angeordnet
  schalter = [0, 0, 0];
    for u = 1:3
        for v = 1:3
            if lambda_vorher(u) == lambda_nachher(v)
                schalter(u) = v;
            end
        end
    end
    Wcopy = W;
    for u = 1:3
        for v = 1:3
            if u == schalter(v)
                W(:,u) = Wcopy(:,v);
            end
        end
    end
%  
%jetzt werden die korrespondierenden Eigenvektoren ausgelesen.
e1= [W(:,1)];
% zweiter eigenwert ist immer 1!
e3= [W(:,3)];

end
  
    


