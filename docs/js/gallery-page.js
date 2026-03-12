document.addEventListener('DOMContentLoaded', () => {
    const galleryGrid = document.getElementById('gallery-grid');
    const modal = document.getElementById('image-modal');
    const modalImg = document.getElementById('modal-image');
    const modalCaption = document.getElementById('modal-caption');
    const modalClose = document.querySelector('.modal-close');

    // Fetch Gallery Data
    fetch('gallery.json')
        .then(response => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            return response.json();
        })
        .then(data => {
            if (galleryGrid) {
                galleryGrid.innerHTML = ''; // Clear loading state

                data.forEach((item, index) => {
                    const card = document.createElement('div');
                    card.className = 'gallery-item animate-on-scroll';
                    card.style.animationDelay = `${index * 0.1}s`; // Staggered animation

                    card.innerHTML = `
                        <div class="gallery-item-img-container">
                            <img src="${item.image}" alt="${item.title}" class="gallery-item-img" loading="lazy">
                        </div>
                        <div class="gallery-item-content">
                            <div class="gallery-item-title">${item.title}</div>
                            ${item.command ? `<div class="gallery-item-command">${item.command}</div>` : ''}
                        </div>
                    `;

                    // Open Modal on Click
                    card.addEventListener('click', () => {
                        modalImg.src = item.image;
                        modalImg.alt = item.title;
                        
                        // Update caption with title and command
                        modalCaption.innerHTML = `
                            <div class="modal-title">${item.title}</div>
                            ${item.command ? `<div class="gallery-item-command">${item.command}</div>` : ''}
                        `;
                        
                        modal.classList.add('active');
                        document.body.style.overflow = 'hidden'; // Prevent scrolling
                    });

                    galleryGrid.appendChild(card);
                });

                // Re-trigger scroll animations
                if (window.IntersectionObserver) {
                    const observer = new IntersectionObserver((entries) => {
                        entries.forEach(entry => {
                            if (entry.isIntersecting) {
                                entry.target.classList.add('visible');
                            }
                        });
                    }, { threshold: 0.1 });
                    
                    document.querySelectorAll('.animate-on-scroll').forEach(el => observer.observe(el));
                }
            }
        })
        .catch(error => {
            console.error('Error loading gallery:', error);
            if (galleryGrid) {
                galleryGrid.innerHTML = '<p style="text-align:center; color:red; grid-column:1/-1;">Failed to load gallery images.</p>';
            }
        });

    // Close Modal Logic
    const closeModal = () => {
        modal.classList.remove('active');
        document.body.style.overflow = ''; // Restore scrolling
        setTimeout(() => {
            modalImg.src = ''; // Clear source to stop memory leak
        }, 300);
    };

    if (modalClose) {
        modalClose.addEventListener('click', closeModal);
    }

    // Close on background click
    if (modal) {
        modal.addEventListener('click', (e) => {
            if (e.target === modal) {
                closeModal();
            }
        });
    }

    // Close on Escape key
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape' && modal.classList.contains('active')) {
            closeModal();
        }
    });
});
