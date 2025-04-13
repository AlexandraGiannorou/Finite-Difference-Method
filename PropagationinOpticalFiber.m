% Δείκτες διάθλασης και παράμετροι
n1 = 1.46; % Δείκτης διάθλασης πυρήνα
n2 = 1.45; % Δείκτης διάθλασης περιβλήματος 
a = 5; % μm, ακτίνα πυρήνα 
lambda = 1550 * 1e-3; % μm
a1 = 0.51; % Παράγοντας σχήματος Crank-Nicolson

% Αναλυτικές τιμές του Neff
neff_analytical = [1.457805787, 1.457161711];

% Πλέγμα
Nx = 150; % Αριθμός κόμβων πλέγματος κατά x
Ny = 150; % Αριθμός κόμβων πλέγματος κατά y
Lx = 4 * a; % μήκος του πλέγματος κατά x
Ly = 4 * a; % μήκος του πλέγματος κατά y
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
dz = 2.5e-6; 
M = Nx * Ny; 
Gamma = 200; % Επαναλήψεις

% Δημιουργία πλέγματος 
x = linspace(-Lx / 2, Lx / 2, Nx);
y = linspace(-Ly / 2, Ly / 2, Ny);
[X, Y] = meshgrid(x, y);

% Δημιουργία πίνακα δείκτη διάθλασης  
ns = n2 * ones(Nx, Ny); % Δίνει στον χώρο τον δείκτη n2
ns(sqrt((X .^ 2) + (Y .^ 2)) <= a) = n1; % Στον χώρο εντός του πυρήνα δίνει τον n1
ns = ns(:);

% Αρχική συνθήκη
E = exp(-((X).^2 + (Y).^2) / (a)^2);
E = E(:); % Μετατροπή σε διάνυσμα

k0 = 2 * pi / lambda;
kref = k0 * n2;
kzero = k0;

% Πενταδιαγώνιος πίνακας S
main_diag_S = (2 * 1i * kref / (1i * dz)) + (2 * a1 / dx^2) + (2 * a1 / dy^2) - a1 * ((kzero^2) * (ns.^2) - (kref^2));
off_diag_x_S = -a1 / dx^2 * ones(M-1, 1);
off_diag_y_S = -a1 / dy^2 * ones(M-Nx, 1);

S = diag(main_diag_S) + diag(off_diag_x_S, 1) + diag(off_diag_x_S, -1) + diag(off_diag_y_S, Nx) + diag(off_diag_y_S, -Nx);
S = sparse(S);

% Πενταδιαγώνιος πίνακας B
main_diag_B = (2 * 1i * kref / (1i * dz)) - (2 * (1-a1) / dx^2) - (2 * (1-a1) / dy^2) + (1-a1) * ((kzero^2) * (ns.^2) - (kref^2));
off_diag_x_B = (1-a1) / dx^2 * ones(M-1, 1);
off_diag_y_B = (1-a1) / dy^2 * ones(M-Nx, 1);

B = diag(main_diag_B) + diag(off_diag_x_B, 1) + diag(off_diag_x_B, -1) + diag(off_diag_y_B, Nx) + diag(off_diag_y_B, -Nx);
B = sparse(B);

% Μηδενισμός στοιχείων πάνω και κάτω από την κύρια διαγώνιο για τον πίνακα S
for i = 1:(Ny-1)
    % Τελευταίο στοιχείο της γραμμής στο πλέγμα
    end_of_row = i * Nx;
    
    if end_of_row < M
        % Μηδενισμός των στοιχείων της πάνω διαγωνίου
        S(end_of_row, end_of_row + 1) = 0;
        
        % Μηδενισμός των στοιχείων της κάτω διαγωνίου
        if end_of_row > 1
            S(end_of_row + 1, end_of_row) = 0;
        end
    end
end

% Μηδενισμός στοιχείων πάνω και κάτω από την κύρια διαγώνιο για τον πίνακα B
for i = 1:(Ny-1)
    % Τελευταίο στοιχείο της γραμμής στο πλέγμα
    end_of_row = i * Nx;
    
    if end_of_row < M
        % Μηδενισμός των στοιχείων της πάνω διαγωνίου
        B(end_of_row, end_of_row + 1) = 0;
        
        % Μηδενισμός των στοιχείων της κάτω διαγωνίου
        if end_of_row > 1
            B(end_of_row + 1, end_of_row) = 0;
        end
    end
end

E_prev = E;

% Διάδοσης σε φανταστική απόσταση
for iter = 1:Gamma/dz
    E_next = S \ (B * E_prev);
    
    % Κανονικοποίηση της δέσμης για σύγκλιση
    E_next = E_next / norm(E_next);

    if iter == Gamma
        break
    end
    E_prev = E_next;
end

En = reshape(E_next, Nx, Ny);
Ep = reshape(E_prev, Nx, Ny);

Eln = log(En ./ Ep) .* (abs(Ep)).^2;
A = abs(Ep);
Elnk = A.^2;

% Ολοκληρώματα
int_E_next = trapz(y, trapz(x, En, 2));
int_E_prev = trapz(y, trapz(x, Ep, 2));
int_ln = trapz(y, trapz(x, Eln));
int_lnk = trapz(y, trapz(x, Elnk));

% Gradient
Egrad = En;
[Egradx, Egrady] = gradient(Egrad, dx, dy); % Υπολογισμός παραγώγων
double_sum_13 = sum(sum((k0^2 * reshape(ns, Nx, Ny).^2 - kref^2) .* abs(Egrad).^2 - (abs(Egradx).^2 + abs(Egrady).^2)));
double_sum_23 = sum(sum(2 * kref * abs(Egrad).^2));

% Υπολογισμός b
b1 = kref + (1 / dz) * log(int_E_next / int_E_prev);
b2 = kref + (1 / dz) * int_ln / int_lnk;
b3 = kref + double_sum_13 / double_sum_23;


% Υπολογισμός ενεργών δεικτών διάθλασης (neff)
neff1 = b1 / k0;
neff2 = b2 / k0;
neff3 = b3 / k0;

% Υπολογισμός σχετικού σφάλματος
relative_error1 = abs((neff1 - neff_analytical) / neff_analytical) * 100;
relative_error2 = abs((neff2 - neff_analytical) / neff_analytical) * 100;
relative_error3 = abs((neff3 - neff_analytical) / neff_analytical) * 100;

% Πλοτάρισμα προφίλ του ρυθμού
figure;
contour(x, y, abs(reshape(E_next, Nx, Ny)).^2, 20); 
colorbar;
title(['|E|^2 at \lambda = ', num2str(lambda * 1e3), ' nm']);
xlabel('x (μm)');
ylabel('y (μm)');

% Αναπαράσταση του ορίου της οπτικής ίνας
hold on;
theta = linspace(0, 2*pi, 100);
x_ina = a * cos(theta);
y_ina = a * sin(theta);
plot(x_ina, y_ina, 'r', 'LineWidth', 1.5); 

disp(['Ενεργοί δείκτες διάθλασης (neff) για το μήκος κύματος ', num2str(lambda * 1e3), ' nm:']);
disp([neff1, neff2, neff3]);
disp(['Σχετικά σφάλματα για το μήκος κύματος ', num2str(lambda * 1e3), ' nm:']);
disp([relative_error1, relative_error2, relative_error3]);
