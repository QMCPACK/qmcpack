
template<typename T, typename Size, typename TT>
shm_ptr_with_raw_ptr_dispatch<T> uninitialized_fill_n(shm_ptr_with_raw_ptr_dispatch<T> first, Size n, TT const& val)
{
  template<class Alloc, typename T, typename Size, typename TT>
  shm_ptr_with_raw_ptr_dispatch<T> uninitialized_fill_n(Alloc & a, shm_ptr_with_raw_ptr_dispatch<T> first, Size n,
                                                        TT const& val)
  {
    template<typename T, typename Size>
    shm_ptr_with_raw_ptr_dispatch<T> destroy_n(shm_ptr_with_raw_ptr_dispatch<T> first, Size n)
    {
      template<class Alloc, typename T, typename Size>
      shm_ptr_with_raw_ptr_dispatch<T> destroy_n(Alloc & a, shm_ptr_with_raw_ptr_dispatch<T> first, Size n)
      {
        template<class It1, typename T, typename Size>
        shm_ptr_with_raw_ptr_dispatch<T> copy_n(It1 first, Size n, shm_ptr_with_raw_ptr_dispatch<T> d_first)
        {
          template<class It1, typename T>
          shm_ptr_with_raw_ptr_dispatch<T> copy(It1 first, It1 last, shm_ptr_with_raw_ptr_dispatch<T> d_first)
          {
            template<class It1, class Size, typename T>
            shm_ptr_with_raw_ptr_dispatch<T> uninitialized_copy_n(It1 f, Size n, shm_ptr_with_raw_ptr_dispatch<T> d)
            {
              template<class Alloc, class It1, class Size, typename T>
              shm_ptr_with_raw_ptr_dispatch<T> uninitialized_copy_n(Alloc & a, It1 f, Size n,
                                                                    shm_ptr_with_raw_ptr_dispatch<T> d)
              {
                template<class It1, typename T>
                shm_ptr_with_raw_ptr_dispatch<T> uninitialized_copy(It1 f, It1 l, shm_ptr_with_raw_ptr_dispatch<T> d);

                template<class T, class Size>
                shm_ptr_with_raw_ptr_dispatch<T> uninitialized_default_construct_n(shm_ptr_with_raw_ptr_dispatch<T> f,
                                                                                   Size n)
                {
                  template<class T, class Size>
                  shm_ptr_with_raw_ptr_dispatch<T> uninitialized_value_construct_n(shm_ptr_with_raw_ptr_dispatch<T> f,
                                                                                   Size n)
                  {
                    template<class Alloc, class T, class Size>
                    shm_ptr_with_raw_ptr_dispatch<T>
                        uninitialized_default_construct_n(Alloc & a, shm_ptr_with_raw_ptr_dispatch<T> f, Size n)
                    {
                      template<class Alloc, class T, class Size>
                      shm_ptr_with_raw_ptr_dispatch<T>
                          uninitialized_value_construct_n(Alloc & a, shm_ptr_with_raw_ptr_dispatch<T> f, Size n)
                      {
                        namespace boost
                        {
                        namespace multi
                        {
                        template<class Alloc, typename T, typename Size>
                        multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> uninitialized_fill_n(
                            Alloc& a,
                            multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> first,
                            Size n,
                            T const& val)
                        {
                          template<class Alloc, typename T, typename Size>
                          multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                              uninitialized_fill(Alloc & a,
                                                 multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                                     first,
                                                 multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                                     last,
                                                 T const& val)
                          {
                            template<class T, class Q1, class Q2, typename Size>
                            multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                copy_n(multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> first,
                                       Size n, multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
                            {
                              template<class T, class ForwardIt, typename Size>
                              multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                  copy_n(ForwardIt first, Size n,
                                         multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
                              {
                                template<class T, class Q1, class Q2, typename Size>
                                multi::array_iterator<T, 1, T*>
                                    copy_n(multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> first,
                                           Size n, multi::array_iterator<T, 1, T*> dest)
                                {
                                  template<class T, class Q1, class Q2>
                                  multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                      copy(multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> first,
                                           multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>> last,
                                           multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
                                  {
                                    template<class T, class ForwardIt>
                                    multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                        copy(ForwardIt first, ForwardIt last,
                                             multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
                                    {
                                      template<class T, class Q1, class Q2>
                                      multi::array_iterator<T, 1, T*>
                                          copy(multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>>
                                                   first,
                                               multi::array_iterator<Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>>
                                                   last,
                                               multi::array_iterator<T, 1, T*> dest)
                                      {
                                        template<class Alloc, class T, class Q, typename Size>
                                        multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                            uninitialized_copy_n(Alloc & a,
                                                                 multi::array_iterator<
                                                                     Q, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q>> first,
                                                                 Size n,
                                                                 multi::array_iterator<
                                                                     T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
                                        {
                                          template<class Alloc, class T, class ForwardIt>
                                          multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                              uninitialized_copy(Alloc & a, ForwardIt first, ForwardIt last,
                                                                 multi::array_iterator<
                                                                     T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>> dest)
                                          {
                                            template<class Alloc, class T, class Q1, class Q2>
                                            multi::array_iterator<T, 1, T*>
                                                uninitialized_copy(Alloc & a,
                                                                   multi::array_iterator<
                                                                       Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>>
                                                                       first,
                                                                   multi::array_iterator<
                                                                       Q1, 1, shm::shm_ptr_with_raw_ptr_dispatch<Q2>>
                                                                       last,
                                                                   multi::array_iterator<T, 1, T*> dest)
                                            {
                                              template<class Alloc, class T, class Size>
                                              multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                                  uninitialized_default_construct_n(
                                                      Alloc & a,
                                                      multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                                          f,
                                                      Size n)
                                              {
                                                template<class Alloc, class T, class Size>
                                                multi::array_iterator<T, 1, shm::shm_ptr_with_raw_ptr_dispatch<T>>
                                                    uninitialized_value_construct_n(
                                                        Alloc & a,
                                                        multi::array_iterator<T, 1,
                                                                              shm::shm_ptr_with_raw_ptr_dispatch<T>> f,
                                                        Size n)
                                                {} // multi
                                              }    // boost
