% Δείκτες διάθλασης και παράμετροι
n_g = double(1.6); % Δείκτης διάθλασης πυρήνα
n_s = double(1.45); % Δείκτης διάθλασης υποστρώματος 
w = double(2); % μm, πλάτος πυρήνα 
h = double(2); % μμ, ύψος πυρήνα
lambda = double(1550 * 1e-3); % μm
a1 = double(0.51); % Παράγοντας σχήματος Crank-Nicolson

% Θεωρητική τιμή του neff από το σαιτ
neff_theoretical = double(1.533532262);

% Πλέγμα
Nx = double(150); % Αριθμός κόμβων πλέγματος κατά x
Ny = double(150); % Αριθμός κόμβων πλέγματος κατά y
Lx = double(4 * w); % μήκος του πλέγματος κατά x
Ly = double(4 * h); % μήκος του πλέγματος κατά y
dx = double(Lx / (Nx - 1));
dy = double(Ly / (Ny - 1));
dz = double(0.2e-6); 
M = double(Nx * Ny); 
Gamma = double(500); 

% Δημιουργία πλέγματος 
x = linspace(-Lx / 2, Lx / 2, Nx);
y = linspace(-Ly / 2, Ly / 2, Ny);
[X, Y] = meshgrid(x, y);

% Δημιουργία πίνακα δείκτη διάθλασης  
ns = n_s * ones(Nx, Ny, 'double'); % Δίνει στον χώρο τον δείκτη n_s
core_region = (abs(X) <= w/2) & (abs(Y) <= h/2);
interface_region_x = (abs(X) == w/2) & (abs(Y) <= h/2);
interface_region_y = (abs(X) <= w/2) & (abs(Y) == h/2);

ns(core_region) = n_g; % Στον χώρο εντός του πυρήνα δίνει τον n_g
ns(interface_region_x | interface_region_y) = (n_g + n_s) / 2; % Μέσος δείκτης διάθλασης στις διαχωριστικές επιφάνειες
ns = ns(:);

% Αρχική συνθήκη
E = exp(-((X).^2 + (Y).^2) / (w)^2);
E = E(:); % Μετατροπή σε διάνυσμα

k0 = double(2 * pi / lambda);
kref = k0 * n_s;
kzero = k0;

% Πενταδιαγώνιος πίνακας S
main_diag_S = (2 * 1i * kref / (1i * dz)) + (2 * a1 / dx^2) + (2 * a1 / dy^2) - a1 * ((kzero^2) * (ns.^2) - (kref^2));
off_diag_x_S = -a1 / dx^2 * ones(M-1, 1, 'double');
off_diag_y_S = -a1 / dy^2 * ones(M-Nx, 1, 'double');

S = diag(main_diag_S) + diag(off_diag_x_S, 1) + diag(off_diag_x_S, -1) + diag(off_diag_y_S, Nx) + diag(off_diag_y_S, -Nx);
S = sparse(S);

% Πενταδιαγώνιος πίνακας B
main_diag_B = (2 * 1i * kref / (1i * dz)) - (2 * (1-a1) / dx^2) - (2 * (1-a1) / dy^2) + (1-a1) * ((kzero^2) * (ns.^2) - (kref^2));
off_diag_x_B = (1-a1) / dx^2 * ones(M-1, 1, 'double');
off_diag_y_B = (1-a1) / dy^2 * ones(M-Nx, 1, 'double');

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
int_E_next = trapz(double(y), trapz(double(x), En, 2));
int_E_prev = trapz(double(y), trapz(double(x), Ep, 2));
int_ln = trapz(double(y), trapz(double(x), Eln));
int_lnk = trapz(double(y), trapz(double(x), Elnk));

%Gradient
Egrad = En;
[Egradx, Egrady] = gradient(Egrad, dx, dy); 
double_sum_13 = sum(sum((k0^2 * reshape(ns, Nx, Ny).^2 - kref^2) .* abs(Egrad).^2 - (abs(Egradx).^2 + abs(Egrady).^2)));
double_sum_23 = sum(sum(2 * kref * abs(Egrad).^2));

%Υπολογισμός β
b1 = kref + (1 / dz) * log(int_E_next / int_E_prev);
b2 = kref + (1 / dz) * int_ln / int_lnk;
b3 = kref + double_sum_13 / double_sum_23;



% Υπολογισμός ενεργών δεικτών διάθλασης (neff)
neff1 = b1 / k0;
neff2 = b2 / k0;
neff3 = b3 / k0;

% Πλοτάρισμα προφίλ
figure;
contour(double(x), double(y), abs(reshape(E_next, Nx, Ny)).^2, 20); % Χρησιμοποιούμε 20 επίπεδα contour
colorbar;
title(['|E|^2 at \lambda = ', num2str(lambda * 1e3), ' nm']);
xlabel('x (μm)');
ylabel('y (μm)');

% Αναπαράσταση κυματοδηγού
hold on;
plot([-w/2, w/2, w/2, -w/2, -w/2], [-h/2, -h/2, h/2, h/2, -h/2], 'r', 'LineWidth', 1.5);

disp(['Ενεργοί δείκτες διάθλασης (neff) για το μήκος κύματος ', num2str(lambda * 1e3), ' nm:']);
disp([neff1, neff2, neff3]);

% Υπολογισμός και εμφάνιση σχετικού σφάλματος σε σύγκριση με τη θεωρητική τιμή
relative_error1 = abs((neff1 - neff_theoretical) / neff_theoretical) * 100;
relative_error2 = abs((neff2 - neff_theoretical) / neff_theoretical) * 100;
relative_error3 = abs((neff3 - neff_theoretical) / neff_theoretical) * 100;

disp(['Σχετικά σφάλματα σε σύγκριση με τη θεωρητική τιμή για το μήκος κύματος ', num2str(lambda * 1e3), ' nm:']);
disp(['Relative Error neff1: ', num2str(relative_error1), '%']);
disp(['Relative Error neff2: ', num2str(relative_error2), '%']);
disp(['Relative Error neff3: ', num2str(relative_error3), '%']);
