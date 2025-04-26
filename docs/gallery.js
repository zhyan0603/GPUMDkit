document.addEventListener('DOMContentLoaded', function() {
    const galleryGrid = document.querySelector('.gallery-grid');
    const prevButton = document.createElement('button');
    const nextButton = document.createElement('button');

    prevButton.className = 'gallery-arrow prev';
    nextButton.className = 'gallery-arrow next';
    prevButton.innerHTML = '❮';
    nextButton.innerHTML = '❯';

    const container = document.createElement('div');
    container.className = 'gallery-grid-container';
    galleryGrid.parentNode.insertBefore(container, galleryGrid);
    container.appendChild(galleryGrid);
    container.appendChild(prevButton);
    container.appendChild(nextButton);

    const scrollAmount = 400;

    prevButton.addEventListener('click', () => {
        galleryGrid.scrollBy({
            left: -scrollAmount,
            behavior: 'smooth'
        });
    });

    nextButton.addEventListener('click', () => {
        galleryGrid.scrollBy({
            left: scrollAmount,
            behavior: 'smooth'
        });
    });

    // Auto scroll functionality
    let scrollInterval;
    const startAutoScroll = () => {
        scrollInterval = setInterval(() => {
            if (galleryGrid.scrollLeft + galleryGrid.clientWidth >= galleryGrid.scrollWidth) {
                galleryGrid.scrollTo({ left: 0, behavior: 'smooth' });
            } else {
                galleryGrid.scrollBy({ left: scrollAmount, behavior: 'smooth' });
            }
        }, 5000);
    };

    const stopAutoScroll = () => {
        clearInterval(scrollInterval);
    };

    // Start auto-scroll by default
    startAutoScroll();

    // Stop auto-scroll when user interacts with the gallery
    container.addEventListener('mouseenter', stopAutoScroll);
    container.addEventListener('mouseleave', startAutoScroll);
}));