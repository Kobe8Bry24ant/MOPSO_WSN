function PlotCosts(pop, rep)

    %pop_costs = [pop.Cost];
    pop_costs = reshape([pop.Cost], 2, []);
    plot(pop_costs(1, :), pop_costs(2, :), 'ko');
    hold on;
    
    %rep_costs = [rep.Cost];
    rep_costs = reshape([rep.Cost], 2, []);
    plot(rep_costs(1, :), rep_costs(2, :), 'r*');
    
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    
    grid on;
    
    hold off;

end