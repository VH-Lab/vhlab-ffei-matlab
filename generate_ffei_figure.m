function [] = generate_ffei_figure(figure_num)
%GENERATE_FFEI_FIGURE: Runs necessary script to generate each figure, given
%the figure number assigned in publication (figure_num).
switch figure_num
    case 1
        triad_ffei;
    case 2
        triad_vs_injection;
    case 3
        triad_Iscaling;
    case 4
        triad_delay;
    case 5
        triad_multilayer;
    case 6
        triad_XOR;
    case 7
        triad_ampagabanmda;
    case 8
        cortex_ffei;
    otherwise
        error('Error: Not a valid figure number. Please enter a figure number ranging from 1-8.')
end
    
end

