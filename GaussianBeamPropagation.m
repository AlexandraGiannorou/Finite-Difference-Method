

lambda0 = 1e-6; 
k0 = 2 * pi / lambda0; 
w0 = 1e-6;
Lx = 12e-6; % υπολογιστικό παράθυρο
Ly = 12e-6; 
Nx = 150; 
Ny = 150; 
dx = Lx / (Nx - 1); % διαστήματα πλέγματος 
dy = Ly / (Ny - 1); 
dz = 1e-6; 
Nz = round(12e-6 / dz); 
a = 0.51; % Crank-Nicolson
relative_threshold = 0.15; % Κάτω όριο τιμής πεδίου
size = (Nx-2) * (Ny-2); % Νέα διάσταση για εσωτερικά σημεία

% Δημιουργία πλέγματος
x = linspace(-Lx / 2, Lx / 2, Nx);
y = linspace(-Ly / 2, Ly / 2, Ny);
[X, Y] = meshgrid(x, y);

% Αρχικοποίηση Γκαουσιανού παλμού
E0 = double(exp(-(X .^ 2 + Y .^ 2) / w0 ^ 2));

% Εξαιρέστε τα άκρα
E = double(E0(2:Nx-1, 2:Ny-1));
E_new = double(zeros(Nx-2, Ny-2));

kref = double(k0);
kzero = double(k0);
n = double(1);

% Πενταδιαγώνιος πίνακας S
S = double(( (2 * 1i * kref / dz) + (2 * a / dx ^ 2) + (2 * a / dy ^ 2) - a * ( (kzero ^ 2) * (n ^ 2) - (kref ^ 2) ) ) * diag(ones(size, 1), 0) ... % Κύρια Διαγώνιος
    - (a / dx ^ 2) * diag(ones(size - 1, 1), 1) - (a / dx ^ 2) * diag(ones(size - 1, 1), -1) ... % Πάνω και κάτω διαγώνιοι
    - (a / dy ^ 2) * diag(ones(size - (Nx-2), 1), Nx-2) - (a / dy ^ 2) * diag(ones(size - (Nx-2), 1), -(Nx-2))); % Απομακρυσμένες
S = sparse(S);

% Πενταδιαγώνιος πίνακας B
B = double(( (2 * 1i * kref / dz) - (2 * (1 - a) / dx ^ 2) - (2 * (1-a) / dy ^ 2) + (1 - a) * ( (kzero ^ 2) * (n ^ 2) - (kref ^ 2) ) ) * diag(ones(size, 1), 0) ...
    + ((1 - a) / dx ^ 2) * diag(ones(size - 1, 1), 1) + ((1 - a) / dx ^ 2) * diag(ones(size - 1, 1), -1) ...
    + ((1 - a) / dy ^ 2) * diag(ones(size - (Nx-2), 1), Nx-2) + ((1 - a) / dy ^ 2) * diag(ones(size - (Nx-2), 1), -(Nx-2)));
B = sparse(B);

% Μηδενισμός στοιχείων πάνω και κάτω από την κύρια διαγώνιο για τον πίνακα S
for i = 1:(Ny-3)
    % Τελευταίο στοιχείο της γραμμής στο πλέγμα
    end_of_row = i * (Nx-2);
    
    if end_of_row < size
        % Μηδενισμός των στοιχείων της πάνω διαγωνίου
        S(end_of_row, end_of_row + 1) = 0;
        
        % Μηδενισμός των στοιχείων της κάτω διαγωνίου
        if end_of_row > 1
            S(end_of_row + 1, end_of_row) = 0;
        end
    end
end

% Μηδενισμός στοιχείων πάνω και κάτω από την κύρια διαγώνιο για τον πίνακα B
for i = 1:(Ny-3)
    % Τελευταίο στοιχείο της γραμμής στο πλέγμα
    end_of_row = i * (Nx-2);
    
    if end_of_row < size
        % Μηδενισμός των στοιχείων της πάνω διαγωνίου
        B(end_of_row, end_of_row + 1) = 0;
        
        % Μηδενισμός των στοιχείων της κάτω διαγωνίου
        if end_of_row > 1
            B(end_of_row + 1, end_of_row) = 0;
        end
    end
end

% Διάδοση του κύματος και αποθήκευση στα ζητούμενα επίπεδα
z_positions = double([1, 2, 5, 10, 12] * 1e-6); % m
E_results = cell(length(z_positions), 1);
z_index = 1;

for n = 1:Nz
    E_vec = reshape(E, size, 1);
    E_new_vec = S\(B * E_vec);
    E_new = reshape(E_new_vec, Nx-2, Ny-2);
    
    E = E_new;

    if abs(n * dz - z_positions(z_index)) < dz / 2
        E_results{z_index} = E;
        z_index = z_index + 1;
    end

    if z_index > length(z_positions)
        break;
    end
end

% Αναλυτική λύση
E_analytical = cell(length(z_positions), 1);
for zi = 1:length(z_positions)
    z = z_positions(zi);
    G = zeros(Nx-2, Ny-2, 'double'); % Εξασφαλίζουμε ότι είναι διπλής ακρίβειας
    z0 = pi * w0^2 / lambda0;
    wz = w0 * sqrt(1 + (lambda0^2 * z^2) / (pi^2 * w0^4));
    Rz = z * (1 + (pi^2 * w0^4) / (z^2 * lambda0^2));

    for xi = 1:(Nx-2)
        for yi = 1:(Ny-2)
            x_val = (xi - (Nx-1)/2) * dx; % Σωστό κεντράρισμα
            y_val = (yi - (Ny-1)/2) * dy; % Σωστό κεντράρισμα
            r = sqrt(x_val^2 + y_val^2);
            G(xi, yi) = (w0 / wz) * exp(-r^2 / wz^2) ...
                * exp(-1i * k0 * r^2 / (2 * Rz)) ...
                * exp(-1i * (k0 * z - atan(z / z0)));
        end
    end
    
    E_analytical{zi} = G;
end

% Ποτάρισμα
for k = 1:length(z_positions)
    % Υπολογιστική Λύση
    figure;
    contour(x(2:Nx-1), y(2:Ny-1), abs(E_results{k}));  % Αναπαράσταση χωρίς τα άκρα
    title(sprintf('Computational solution at z = %.1f μm', z_positions(k) * 1e6));
    xlabel('x (m)');
    ylabel('y (m)');
    colorbar;

    % Αναλυτική Λύση
    figure;
    contour(x(2:Nx-1), y(2:Ny-1), abs(E_analytical{k}));  % Αναπαράσταση χωρίς τα άκρα
    title(sprintf('Analytical solution at z = %.1f μm', z_positions(k) * 1e6));
    xlabel('x (m)');
    ylabel('y (m)');
    colorbar;
    
    % Υπολογισμός ορίου βάση της μέγιστης έντασης του πεδίου
    max_intensity = max(abs(E_results{k}(:)));
    dynamic_threshold = relative_threshold * max_intensity;
    
    % Υπολογισμός ακρίβειας όπου το πεδίο είναι μεγαλύτερο απο το όριο 
    overv = abs(E_results{k}) > dynamic_threshold;
    Tv = abs(E_analytical{k}(overv));
    Ov = abs(E_results{k}(overv));
    accuracy = mean(100 - ((abs(Tv - Ov) ./ Tv) * 100));
    fprintf('Mean accuracy at z = %.1f μm (above dynamic threshold): %.2f%%\n', z_positions(k) * 1e6, accuracy);
end
